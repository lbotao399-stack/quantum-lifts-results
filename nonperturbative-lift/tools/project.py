#!/usr/bin/env python3
"""Navigate, verify and replay the local non perturbative lift project (stdlib)."""
import argparse
import collections
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from urllib.parse import unquote, urlsplit
import zipfile

ROOT = Path(__file__).resolve().parents[1]
ZIP = ROOT / 'artifacts/nonperturbative-certificates.zip'
ARCHIVE_PIN = '2d94d52cbd29e7f498630158ad92b5464992cc7878002860fee68791f478b066'
ANCHOR = re.compile(r'<a id="([\w.-]+)"></a>')
STATES = ('proved', 'conditional', 'literature-input', 'finite-evidence', 'conjecture', 'open')


def load(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8'))


def require(test, message):
    if not test:
        raise ValueError(message)


def digest(stream):
    h = hashlib.sha256()
    while chunk := stream.read(1024*1024):
        h.update(chunk)
    return h.hexdigest()


def hash_file(path):
    with path.open('rb') as f:
        return digest(f)


def inside(base, name):
    path = (base / name).resolve()
    require(not Path(name).is_absolute() and path.is_relative_to(ROOT), f'Path outside project: {name}')
    return path


def artifact_index():
    a = load('artifacts/inventory.json')
    require(a['archive_sha256'] == ARCHIVE_PIN and hash_file(ZIP) == ARCHIVE_PIN, 'Archive SHA-256 mismatch')
    rows = {r['path']: r for r in a['payloads']}
    require(len(rows) == len(a['payloads']), 'Duplicate inventory member')
    return a, rows


def member_bytes(z, name, rows):
    require(name in rows, f'Unknown certificate payload: {name}')
    data = z.read(name)  # zipfile checks the CRC as it reads to EOF.
    ref = rows[name]
    require(len(data) == ref['bytes'] and hashlib.sha256(data).hexdigest() == ref['sha256'], f'Bad payload: {name}')
    return data


def dump(data, path=None):
    text = json.dumps(data, ensure_ascii=False, indent=2) + '\n'
    if path:
        out = Path(path)
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=out.parent, delete=False) as f:
            f.write(text)
        os.replace(f.name, out)
    else:
        print(text, end='')


def continuation(c):
    key = c['id']
    path = 'OPEN_PROBLEMS.md'
    require(c.get('continuation') == f'{path}#{key}', f'No continuation indexed for {key}. Use claims --status open.')
    text = (ROOT / path).read_text(encoding='utf-8')
    marker = f'<a id="{key}"></a>'
    require(text.count(marker) == 1, f'Missing or ambiguous continuation anchor: {key}')
    block = text.split(marker, 1)[1].lstrip()
    # Each task has a heading. A group heading can follow the last task in a group.
    heading = re.match(r'^(#{2,3}) [^\n]+\n', block)
    require(heading is not None, f'Missing task heading: {key}')
    end = re.search(r'^(?:<a id="|#{1,3} )', block[heading.end():], flags=re.M)
    return block[:heading.end() + end.start()].strip() if end else block.strip()


def print_claim(c):
    print(f"{c['id']} [{c['status']}]  {c['title']}")
    for key in ('theory','sector','coefficient_ring','hypotheses','statement','derivation','evidence','replays','depends_on','limits'):
        v = c[key]
        print(f"  {key}: " + (', '.join(v) or 'none' if isinstance(v,list) else v))
    if c.get('continuation'): print(f"  continuation: {c['continuation']} (show {c['id']} --next)")


def replay(args):
    catalog = load('checks/catalog.json')
    if args.list:
        for key,s in catalog.items(): print(f"{key:22} {s['requires']:7} {s['scope']}")
        return 0
    ids = list(catalog) if args.id == 'all' else [args.id]
    require(all(i in catalog for i in ids), 'Unknown replay. Use replay --list.')
    _, rows = artifact_index()
    originals = load('sources/inventory.json')['files']
    for name,rec in originals.items():
        require(hash_file(inside(ROOT,name))==rec['sha256'], f'Upstream file changed: {name}')
    env = {k:v for k,v in os.environ.items() if k not in ('PYTHONOPTIMIZE','PYTHONPATH')}
    result = {'timestamp_utc':dt.datetime.now(dt.timezone.utc).isoformat(), 'python':platform.python_version(), 'implementation':platform.python_implementation(), 'assertions_enabled_in_children':True, 'archive_sha256':ARCHIVE_PIN, 'results':[]}
    import importlib.metadata
    for package in ('sympy','numpy'):
        try: result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError: result[package] = None
    for key in ids:
        s = catalog[key]; kind = s['kind']; name = s['path']; start = time.monotonic()
        print(f'Running {key} ...', flush=True)
        with tempfile.TemporaryDirectory(prefix='np-lift-') as tmp:
            wd = Path(tmp)
            for d in s.get('mkdir',[]): (wd/d).mkdir(parents=True, exist_ok=True)
            if kind == 'archived':
                target = wd/Path(name).name
                with zipfile.ZipFile(ZIP) as z: target.write_bytes(member_bytes(z, name, rows))
            else:
                src = inside(ROOT,name); target = src
                if kind == 'original':
                    target = wd/src.name; shutil.copyfile(src,target)
                    for inp in s.get('inputs',[]): shutil.copyfile(inside(ROOT,inp), wd/Path(inp).name)
            if kind == 'function':
                cmd = ['-c', 'import sys; sys.path.insert(0,sys.argv[1]); import check; getattr(check,sys.argv[2])()',str(target.parent),s['function']]
            else:
                cmd = [str(target)]
            proc = subprocess.run([sys.executable,'-B',*cmd],cwd=wd,env=env,capture_output=True,text=True,errors='replace')
            def clean(text): return text.strip().replace(str(wd),'<replay-dir>').replace(str(ROOT),'<project>')
            row = {'id':key,'status':'passed' if proc.returncode==0 else 'failed','returncode':proc.returncode,'scope':s['scope'],'script_sha256':hash_file(target),'seconds':round(time.monotonic()-start,3),'stdout':clean(proc.stdout),'stderr':clean(proc.stderr)}
            print(f"{row['status']} ({row['seconds']}s)" + (f": {row['stderr']}" if row['returncode'] else ''),flush=True)
            result['results'].append(row)
    result['passed'] = sum(r['status']=='passed' for r in result['results']); result['total'] = len(ids)
    if args.report: dump(result,args.report)
    elif len(ids)==1: print(result['results'][0]['stdout'])
    print(f"{result['passed']}/{result['total']} passed")
    return int(result['passed'] != len(ids))


def markdown_links(text):
    # In these documents Markdown links occur in prose/tables. TeX has legitimate ]( constructs.
    math = code = False
    for line in text.splitlines():
        body = line.removeprefix('>').strip()
        if body.startswith('```'): code = not code; continue
        if code: continue
        if body == '$$': math = not math; continue
        if math: continue
        body = re.sub(r'`[^`]*`|(?<!\\)\$(?:\\.|[^$])*?\$', '', body)
        for m in re.finditer(r'\[[^\]\n]*\]\((?:<([^>]+)>|([^\s)]+))(?:\s+"[^"]*")?\)', body):
            yield m.group(1) or m.group(2)


def verify_publication(excluded, asset_dir, manuscript_source_sha256):
    publication = load('artifacts/publication.json')
    witnesses = {w['path']: w for w in publication['witnesses']}
    require(len(witnesses) == len(excluded) == len(publication['witnesses']), 'Witness count mismatch')
    checked = 0
    for record in excluded:
        name = record['name']
        witness = witnesses[name]
        require(witness['bytes'] == record['bytes'], f'Historical size mismatch: {name}')
        if 'sha256_from_verified_status' in record:
            require(witness['sha256'] == record['sha256_from_verified_status'], f'Historical digest mismatch: {name}')
        if 'file' in witness:
            path = inside(ROOT, witness['file'])
        elif asset_dir is not None:
            path = asset_dir / witness['release_asset']
        else:
            continue
        require(path.stat().st_size == witness['gzip_bytes'] and hash_file(path) == witness['gzip_sha256'], f'Bad gzip file: {name}')
        with gzip.open(path, 'rb') as stream:
            require(digest(stream) == witness['sha256'] and stream.tell() == witness['bytes'], f'Bad polynomial bytes: {name}')
        checked += 1

    history = publication['historical_calculations']
    index = load(history['index'])
    files = {r['path']: r for r in index['files']}
    require(len(files) == len(index['files']) == history['members'], 'Historical file count mismatch')
    manuscript = index['pre_split_manuscript']
    require(manuscript['historical_source_sha256'] == manuscript_source_sha256, 'Historical manuscript provenance mismatch')
    require(files[manuscript['path']]['sha256'] == manuscript['public_text_sha256'], 'Published manuscript digest mismatch')
    history_checked = 0
    if asset_dir is not None:
        path = asset_dir / history['release_asset']
        require(path.stat().st_size == history['bytes'] and hash_file(path) == history['sha256'], 'Historical ZIP digest mismatch')
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            require(len(names) == len(set(names)) and set(names) == set(files) | {'manifest.json'}, 'Historical ZIP members differ')
            require(json.loads(archive.read('manifest.json')) == index, 'Historical ZIP manifest differs')
            for name, record in files.items():
                require(archive.getinfo(name).file_size == record['bytes'], f'Historical file size mismatch: {name}')
                with archive.open(name) as member:
                    require(digest(member) == record['sha256'], f'Historical file digest mismatch: {name}')
                history_checked += 1
    return {'published_polynomials_checked': checked, 'historical_files_indexed': len(files),
            'historical_files_checked': history_checked}


def verify(args):
    counts = {}; sources = load('sources/inventory.json'); files = sources['files']; anchors = {}
    for name,r in files.items():
        p=inside(ROOT,name);require(p.is_file() and p.stat().st_size == r['bytes'] and hash_file(p)==r['sha256'],f'Source mismatch: {name}')
    counts['pinned_sources'] = len(files)
    a,rows = artifact_index(); counts['archive_payloads'] = len(rows)
    with zipfile.ZipFile(ZIP) as z:
        names=z.namelist();require(len(names)==len(set(names)) and set(names)==set(rows)|{'manifest.json','README.md'},'ZIP index differs')
        m=json.loads(z.read('manifest.json'));z.read('README.md')
        require(len(m['included'])==len(rows) and {r['path']:{k:r[k] for k in ('path','bytes','sha256')} for r in m['included']}==rows,'Internal/external manifests differ')
        excluded=[{**{k:v for k,v in r.items() if k!='source'},'name':r['source'].split('/work/',1)[-1]} for r in m['archived_local']]
        require(excluded==a['not_in_archive'] and m['missing']==a['recorded_missing'], 'Excluded records differ')
        expected = set()
        for name in rows:
            require(not PurePosixPath(name).is_absolute() and '..' not in PurePosixPath(name).parts,'Unsafe ZIP member')
            member_bytes(z,name,rows)
            if name.endswith(('.json','.md')):
                p=ROOT/'evidence'/name;expected.add(p);require(p.is_file() and hash_file(p)==rows[name]['sha256'],f'Expanded evidence mismatch: {name}')
        require(expected=={p for p in (ROOT/'evidence').rglob('*') if p.is_file()},'Extra or missing evidence file')
    counts.update(expanded_evidence=len(expected),witnesses_outside_certificate_zip=len(excluded))
    counts.update(verify_publication(excluded, args.release_dir, sources['manuscript_source_sha256']))
    def ref(base,link):
        if link.startswith('archive:'):
            require(link[8:] in rows,'Missing '+link);return
        url=urlsplit(link)
        require(not url.scheme and not url.netloc, 'Expected local reference: '+link)
        p=inside(base,unquote(url.path or '.'));require(p.exists(), 'Missing reference: '+str(base.relative_to(ROOT))+'/'+link)
        if url.fragment:
            require(p.is_file(),f'Fragment on non-file: {link}')
            if p not in anchors: anchors[p]=set(ANCHOR.findall(p.read_text()))
            require(unquote(url.fragment) in anchors[p],f'Missing explicit anchor: {link}')
    sect=load('research/sections.json'); secs=sect['sections']; require(sect['source_sha256']==sources['manuscript_source_sha256'], 'Manuscript provenance mismatch')
    require([s['number'] for s in secs]==list(range(27)), 'Chapter order changed')
    pinned_sections = set(); eqs = blocks = tables = 0
    for s in secs:
        p=inside(ROOT,s['file']);require(hash_file(p)==s['sha256'],f'Section digest mismatch: {p.name}')
        t=p.read_text(encoding='utf-8');lines=t.splitlines(); ids=ANCHOR.findall(t)
        require(len(ids)==len(set(ids)), f'Duplicate anchor: {p.name}')
        headings=[(len(m[1]),m[2]) for l in lines if (m:=re.fullmatch(r'(#{1,6}) (.+)',l))]
        require(headings==[(h['level'],h['title']) for h in s['headings']],f'Changed headings: {p.name}')
        for h in s['headings']:
            require('<a id="'+h['id']+'"></a>\n\n'+'#'*h['level']+' '+h['title']+'\n' in t,'Unpaired heading: '+h['id'])
        delim=sum(l.lstrip('> ').strip()=='$$' for l in lines)
        nblocks=sum(l.startswith('> **') and (i==0 or not lines[i-1].startswith('>')) for i,l in enumerate(lines))
        ntables=sum(l.startswith('|') and (i==0 or not lines[i-1].startswith('|')) for i,l in enumerate(lines))
        require(delim==2*s['equations'] and nblocks==s['proof_blocks'] and ntables==s['tables'],f'Section statistics mismatch: {p.name}')
        pinned_sections.add(p);eqs+=delim//2;blocks+=nblocks;tables+=ntables
    require(pinned_sections == set((ROOT/'research').glob('[0-9][0-9]-*.md')), 'Unindexed section')
    counts.update(numbered_sections=len(secs)-1,display_equations=eqs,proof_or_example_blocks=blocks,tables=tables)
    catalog=load('checks/catalog.json')
    for k,v in catalog.items():
        require(v['kind'] in ('archived','project','original','function') and v['scope'] and v['requires'] in ('stdlib','sympy','sympy+numpy'),'Bad replay: '+k)
        ref(ROOT,('archive:' if v['kind']=='archived' else '')+v['path'])
        for x in v.get('inputs',[]): ref(ROOT,x)
    data=load('claims.json'); claims=data['claims']; keys={c['id']:c for c in claims};require(len(keys)==len(claims) and set(STATES)==set(data['statuses']),'Invalid claim index or statuses')
    visited=set()
    def walk(k,path):
        require(k in keys, 'Missing dependency: '+k);require(k not in path,'Dependency cycle at '+k)
        if k not in visited:
            for dep in keys[k]['depends_on']:walk(dep,path|{k})
            visited.add(k)
    statuses=collections.Counter()
    for c in claims:
        key=c['id'];require(re.fullmatch(r'[a-z0-9]+(?:-[a-z0-9]+)*',key) and c['status'] in STATES,'Bad claim: '+key)
        for field in ('title','statement','theory','sector','coefficient_ring','hypotheses','limits','derivation'):require(isinstance(c.get(field),str) and c[field].strip(),f'{key}: missing {field}')
        for field in ('evidence','replays','depends_on'):require(isinstance(c.get(field),list),f'{key}: invalid {field}')
        require(c['evidence'] and '#' in c['derivation'],f'{key}: missing primary record/anchor');ref(ROOT,c['derivation'])
        for r in c['evidence']:ref(ROOT,r)
        for r in c['replays']:require(r in catalog,f'{key}: bad replay {r}')
        if c['status']=='open' or c.get('continuation'):
            continuation(c);ref(ROOT,c['continuation'])
        walk(key,set());statuses[c['status']]+=1
    docs={p for p in ROOT.rglob('*.md') if str(p.relative_to(ROOT)) not in files and not p.is_relative_to(ROOT/'evidence')}
    links=0
    for p in docs:
        for url in markdown_links(p.read_text()):
            if urlsplit(url).scheme in ('http','https'):continue
            ref(p.parent,(p.name if url.startswith('#') else '')+url);links+=1
    counts.update(claims=len(claims),claims_by_status=dict(sorted(statuses.items())),registered_replays=len(catalog),authored_documents=len(docs),checked_local_markdown_links=links)
    report={'status':'passed','checks':'Source and chapter digests, certificate records, claim references, local links, and published polynomial bytes. With release files supplied, the calculation archive and every indexed member are also checked.','archive_sha256':ARCHIVE_PIN,'counts':counts}
    if args.report: dump(report,args.report)
    dump(report);return 0


def main():
    pa=argparse.ArgumentParser(description=__doc__); sub=pa.add_subparsers(dest='cmd',required=True)
    a=sub.add_parser('claims',help='Search the claim ledger; all search words must match');a.add_argument('query',nargs='*');a.add_argument('--status',choices=STATES);a.add_argument('--json',action='store_true')
    a=sub.add_parser('show',help='Open one claim, including limits and evidence');a.add_argument('id');fmt=a.add_mutually_exclusive_group();fmt.add_argument('--json',action='store_true');fmt.add_argument('--next',action='store_true',help='Read only this task from OPEN_PROBLEMS.md')
    a=sub.add_parser('verify',help='Check project integrity (no optional packages)');a.add_argument('--report')
    a.add_argument('--release-dir', type=Path, help='Also check both downloaded release assets and their uncompressed contents')
    a=sub.add_parser('extract',help='Copy one hash-verified archived payload to a new path');a.add_argument('name');a.add_argument('output')
    a=sub.add_parser('replay',help='Run calculations with enabled assertions in isolated temporary directories');a.add_argument('id',nargs='?');a.add_argument('--list',action='store_true');a.add_argument('--report')
    args=pa.parse_args()
    if args.cmd=='replay':return replay(args)
    if args.cmd=='verify':return verify(args)
    if args.cmd=='extract':
        _,rows=artifact_index()
        with zipfile.ZipFile(ZIP) as z:data=member_bytes(z,args.name,rows)
        with Path(args.output).open('xb') as f:f.write(data)
        print(f"{rows[args.name]['sha256']}  {args.output}");return 0
    claims=load('claims.json')['claims']
    if args.cmd=='show':
        matches=[c for c in claims if c['id']==args.id];require(len(matches)==1,'Unknown claim: '+args.id)
        if args.next: print(continuation(matches[0]))
        elif args.json: dump(matches[0])
        else: print_claim(matches[0])
        return 0
    def search_text(c):
        text=' '.join(str(v) for v in c.values()).casefold()
        return text+' '+re.sub('[(){}_-]','',text)
    matches=[c for c in claims if (not args.status or args.status==c['status']) and all(word.casefold() in search_text(c) for word in args.query)]
    if args.json:dump(matches)
    else:
        for c in matches: print(f"{c['id']:32} {c['status']:17} {c['title']}")
        print(f'{len(matches)} claim(s)')
    return 0


if __name__=='__main__':
    try:sys.exit(main())
    except (OSError,ValueError,KeyError,zipfile.BadZipFile) as exc:sys.exit(f'Error: {exc}')
