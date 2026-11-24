#!/usr/bin/env python3
"""Freeze the public claim-source graph without importing private experiments.

Only tracked public documentation at an explicit immutable revision, explicit
links from that documentation, and the fetched issue snapshot are considered.
This inventories statements; it does not promote them to established results.
"""
import argparse
import gzip
import hashlib
import json
import posixpath
import re
import subprocess
from collections import deque
from pathlib import Path
from urllib.parse import unquote, urlsplit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASELINE = '19cd8718b752b04104d6ff0becee42d580b27ca3'
REPOSITORY = 'lbotao399-stack/quantum-lifts-results'
BLOCKED = ('/isolated/', 'highly_innovative', '/quarantine/', '/private/')


def command(*args, binary=False):
    try:
        result = subprocess.run(['git', '-C', str(ROOT), *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30)
    except subprocess.TimeoutExpired as error:
        raise ValueError('Git object lookup timed out; reference is unverified: ' + ' '.join(args)) from error
    if result.returncode:
        raise ValueError(result.stderr.decode(errors='replace').strip())
    return result.stdout if binary else result.stdout.decode()


def permitted(path):
    return not any(part in ('/' + path).lower() for part in BLOCKED)


def doc(path):
    return path.endswith(('.md', '.tex')) and Path(path).name != 'AGENTS.md' and '/runtime/' not in path


def text_links(body):
    # A closing bracket followed by a parenthesis is common mathematical
    # notation. Require the complete Markdown link, then a path-like target.
    values = re.findall(r'(?<!!)(?<!\\)\[[^\[\]\n]+\]\(([^\s()]+)(?:\s+"[^"]*")?\)', body)
    values = [value for value in values
              if not any(character in value for character in '\\{}$')
              and (value.startswith(('#', 'https://', 'http://')) or '/' in value or '.' in value)]
    values += re.findall(r'https://[^\s<>"}]+', body)
    result = []
    for value in values:
        value = value.split(' "', 1)[0].strip()
        previous = None
        while previous != value:
            previous = value
            value = value.rstrip('.,;:!?*`\'')
            if value.endswith(')') and value.count(')') > value.count('('):
                value = value[:-1]
        if value and not value.startswith(('mailto:', 'data:')):
            result.append(value)
    return sorted(set(result))


def resolve(url, commit, path):
    if url.startswith('#'):
        return None
    parsed = urlsplit(url)
    if parsed.scheme:
        match = re.match(r'^/' + re.escape(REPOSITORY) + r'/(blob|tree)/([^/]+)/(.*)', parsed.path)
        if parsed.netloc == 'github.com' and match:
            kind, revision, item = match.groups()
            if revision in ('main', 'HEAD'):
                revision = BASELINE
            if not re.fullmatch('[0-9a-f]{40}', revision):
                return {'kind': 'external', 'url': url, 'reason': 'non-immutable revision requires explicit metadata'}
            return {'kind': 'git', 'commit': revision, 'path': unquote(item), 'anchor': parsed.fragment}
        return {'kind': 'external', 'url': url}
    item = posixpath.normpath(posixpath.join(posixpath.dirname(path), unquote(parsed.path)))
    return {'kind': 'git', 'commit': commit, 'path': item, 'anchor': parsed.fragment}


def write_gzip(path, data):
    raw = json.dumps(data, ensure_ascii=False, sort_keys=True, separators=(',', ':')).encode()
    with path.open('wb') as output:
        with gzip.GzipFile(filename='', fileobj=output, mode='wb', mtime=0) as archive:
            archive.write(raw)


def build():
    paths = command('ls-tree', '-r', '--name-only', BASELINE).splitlines()
    roots = [p for p in paths if doc(p) and permitted(p) and not p.startswith('research-roadmap/')]
    queue = deque((BASELINE, p) for p in roots)
    snapshot = json.loads((HERE / 'issue-1-snapshot.json').read_text())
    sources = {}
    references = {}
    missing = []
    exclusions = []
    objects = {}

    def inspect_reference(reference, origin):
        key = json.dumps(reference, sort_keys=True)
        if key in references:
            references[key]['occurrences'].append(origin)
            return
        references[key] = {**reference, 'occurrences': [origin]}
        if reference['kind'] != 'git':
            return
        commit, path = reference['commit'], reference['path']
        if not permitted(path):
            exclusions.append({'commit': commit, 'path': path, 'reason': 'quarantine exclusion; contents not read'})
            references[key]['state'] = 'excluded'
            return
        try:
            kind = command('cat-file', '-t', commit + ':' + path).strip()
            object_id = command('rev-parse', commit + ':' + path).strip()
        except ValueError as error:
            references[key]['state'] = 'missing'
            missing.append({'reference': reference, 'origin': origin, 'reason': str(error)})
            return
        references[key].update(state='resolved', git_object=object_id, object_type=kind)
        objects[commit + ':' + path] = object_id
        if kind == 'tree':
            # A directory link advertises its public result documentation.
            for child in command('ls-tree', '-r', '--name-only', commit, '--', path).splitlines():
                if doc(child) and permitted(child):
                    queue.append((commit, child))
        elif doc(path):
            queue.append((commit, path))

    for entry in [{'id': 'issue-body', 'body': snapshot['issue']['body']}, *snapshot['comments']]:
        identifier = 'issue-1:' + str(entry['id'])
        sources[identifier] = {'kind': 'issue-snapshot', 'id': identifier, 'text': entry['body'],
                               'sha256': hashlib.sha256(entry['body'].encode()).hexdigest()}
        for url in text_links(entry['body']):
            reference = resolve(url, BASELINE, 'RESEARCH_DISCUSSION.md')
            if reference:
                inspect_reference(reference, identifier)

    seen = set()
    checkpoint = HERE / 'inventory-checkpoint.json.gz'
    if checkpoint.exists():
        saved = json.loads(gzip.decompress(checkpoint.read_bytes()))
        if saved['baseline'] != BASELINE:
            raise ValueError('Inventory checkpoint belongs to a different baseline.')
        sources = saved['sources']; references = saved['references']
        missing = saved['missing']; exclusions = saved['exclusions']; objects = saved['objects']
        seen = set(tuple(x) for x in saved['seen']); queue = deque(tuple(x) for x in saved['queue'])
        print(json.dumps({'resumed_documents': len(seen), 'queued': len(queue)}), flush=True)
    while queue:
        commit, path = queue.popleft()
        if (commit, path) in seen:
            continue
        seen.add((commit, path))
        if not permitted(path):
            continue
        raw = command('show', commit + ':' + path, binary=True)
        body = raw.decode()
        identifier = commit + ':' + path
        sources[identifier] = {'kind': 'git-document', 'id': identifier, 'commit': commit, 'path': path,
                               'git_blob': command('rev-parse', identifier).strip(),
                               'sha256': hashlib.sha256(raw).hexdigest(), 'text': body}
        for url in text_links(body):
            reference = resolve(url, commit, path)
            if reference:
                inspect_reference(reference, identifier)
        if len(seen) % 25 == 0:
            write_gzip(checkpoint, {'baseline': BASELINE, 'sources': sources, 'references': references,
                                   'missing': missing, 'exclusions': exclusions, 'objects': objects,
                                   'seen': sorted(seen), 'queue': list(queue)})
            print(json.dumps({'documents_read': len(seen), 'queued': len(queue), 'missing_references': len(missing)}), flush=True)
    output = HERE / 'inventory'
    output.mkdir(exist_ok=True)
    write_gzip(output / 'public_sources.json.gz', sources)
    manifest = {'schema_version': 1, 'repository': REPOSITORY, 'baseline': BASELINE,
                'scope': 'All public tracked documentation at the baseline and its reachable immutable repository-document links, plus issue 1 snapshot; not the entire history of unrelated branches.',
                'root_documents': roots, 'document_count': len(sources), 'git_document_count': len(seen),
                'reference_count': len(references), 'missing_reference_count': len(missing),
                'issue_comments': len(snapshot['comments']),
                'sources_sha256': hashlib.sha256((output / 'public_sources.json.gz').read_bytes()).hexdigest(),
                'references': list(references.values()), 'missing': missing, 'exclusions': exclusions}
    (output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({key: manifest[key] for key in ('document_count', 'reference_count', 'missing_reference_count', 'issue_comments')}))


if __name__ == '__main__':
    build()
