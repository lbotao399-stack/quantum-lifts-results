#!/usr/bin/env python3
"""Validate and render the versioned research rubric; never certify physics by fiat."""
import argparse
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
GENERATED = ('SCORECARD.md', 'MILESTONES.md', 'EVIDENCE.md')


def validate(data, evidence):
    def require(condition, message):
        if not condition:
            raise ValueError(message)

    require(data['schema_version'] == 1, 'unsupported schema version')
    rows = data['milestones']
    tracks = {t['id']: t for t in data['tracks']}
    by_id = {r['id']: r for r in rows}
    require(len(rows) == data['milestone_count'] == 60, 'expected 60 milestones')
    require(len(by_id) == len(rows), 'duplicate milestone ID')
    require(len(tracks) == len(data['tracks']) == 7, 'duplicate or missing track')
    require(sum(t['maximum'] for t in tracks.values()) == data['total_points'] == 100,
            'track maxima must total 100')
    require(re.fullmatch(r'[0-9a-f]{40}', data['baseline_commit']), 'unpinned baseline')

    for row in rows:
        identifier = row['id']
        require(re.fullmatch(r'[A-G][0-9]{2}', identifier), f'invalid ID: {identifier}')
        require(row['track'] in tracks and identifier[0] == row['track'], f'bad track: {identifier}')
        require(type(row['points']) is int and row['points'] > 0, f'invalid weight: {identifier}')
        require(row['status'] in ('met', 'partial', 'open'), f'invalid status: {identifier}')
        for key in ('title', 'scope', 'current_assessment', 'next_action', 'required_artifact', 'audit_basis'):
            require(isinstance(row[key], str) and row[key].strip(), f'missing {key}: {identifier}')
        require(len(row['acceptance']) >= 2 and all(isinstance(s,str) and s.strip() for s in row['acceptance']),
                f'incomplete acceptance gate: {identifier}')
        require(row['evidence'] and all(e in evidence for e in row['evidence']), f'unknown evidence: {identifier}')
        require(len(set(row['dependencies'])) == len(row['dependencies']), f'duplicate dependency: {identifier}')
        require(all(d in by_id and d != identifier for d in row['dependencies']), f'bad dependency: {identifier}')
        if row['status'] == 'met':
            require(row['audit_basis'] in ('published_certificate_audit', 'fresh_complete_replay'),
                    f'missing met-gate provenance: {identifier}')
            require(all(by_id[d]['status'] == 'met' for d in row['dependencies']),
                    f'met milestone has unmet dependency: {identifier}')

    for key, item in evidence.items():
        require(item['url'].startswith('https://'), f'non-HTTPS evidence: {key}')
        require(item.get('finding') and item.get('review'), f'incomplete evidence record: {key}')
        if 'commit' in item:
            require(re.fullmatch(r'[0-9a-f]{40}', item['commit']), f'unpinned evidence: {key}')
            require('/blob/' + item['commit'] + '/' + item['path'] in item['url'], f'inconsistent evidence path: {key}')
        require(not any(x in item.get('path','').lower() for x in ('/isolated/', 'highly_innovative')),
                f'quarantined evidence must not be published here: {key}')

    done, active = set(), set()
    def visit(identifier):
        require(identifier not in active, f'dependency cycle through {identifier}')
        if identifier in done:
            return
        active.add(identifier)
        for dep in by_id[identifier]['dependencies']:
            visit(dep)
        active.remove(identifier)
        done.add(identifier)
    for identifier in by_id:
        visit(identifier)
    for identifier, track in tracks.items():
        require(sum(r['points'] for r in rows if r['track'] == identifier) == track['maximum'],
                f'track weight mismatch: {identifier}')
    require(sum(r['points'] for r in rows) == 100, 'milestone weights must total 100')
    require(len(set(data['priority_queue'])) == len(data['priority_queue']), 'duplicate queue item')
    require(all(i in by_id and by_id[i]['status'] != 'met' for i in data['priority_queue']), 'invalid next-action queue')
    return by_id


def earned(row):
    return row['points'] if row['status'] == 'met' else 0


def totals(data):
    rows = data['milestones']
    return dict(score=sum(map(earned, rows)), maximum=data['total_points'],
                met=sum(r['status'] == 'met' for r in rows), count=len(rows),
                tracks=[dict(id=t['id'], title=t['title'],
                             earned=sum(earned(r) for r in rows if r['track'] == t['id']),
                             maximum=t['maximum'],
                             met=sum(r['status'] == 'met' for r in rows if r['track'] == t['id']),
                             count=sum(r['track'] == t['id'] for r in rows)) for t in data['tracks']])


def render(data, evidence):
    total = totals(data)
    by_id = {r['id']:r for r in data['milestones']}
    stamp = f"Rubric {data['rubric_version']} · {data['assessment_date']} ({data['timezone']})"
    pin = data['baseline_commit']
    base = f"https://github.com/{data['repository']}/commit/{pin}"
    prelude = ['<!-- Generated by score.py from milestones.json and evidence.json. -->', '']
    card = prelude + ['# Current score', '', f"**{total['score']}/100 points; {total['met']}/{total['count']} milestones met.**", '',
        stamp, '', f'Public evidence baseline: [{pin}]({base}).', '',
        'This is a weighted completion score for the declared research programme. It is not a percentage of the infinite BPS spectrum, a probability of success, or a measure of elapsed work. Only a fully met gate earns its integer weight. Partial and open gates earn zero.', '',
        'The assessment audits published certificates and their stated replay scope. It does not claim that every large calculation was rerun for this assessment. The roadmap itself earns no research points.', '',
        '| Track | Earned | Available | Met milestones |', '|---|---:|---:|---:|']
    card += [f"| {t['id']} · {t['title']} | {t['earned']} | {t['maximum']} | {t['met']}/{t['count']} |" for t in total['tracks']]
    card += [f"| **Total** | **{total['score']}** | **100** | **{total['met']}/{total['count']}** |", '',
             '## Exact score calculation', '',
             '`' + ' + '.join(str(t['earned']) for t in total['tracks']) + f" = {total['score']}`", '',
             'Earned gates: ' + ', '.join(f"[{r['id']}](MILESTONES.md#{r['id'].lower()}) ({r['points']})" for r in data['milestones'] if r['status']=='met') + '.', '',
             '## All milestones', '',
             '| ID | Milestone | Earned / maximum | Status |', '|---|---|---:|---|']
    card += [f"| [{r['id']}](MILESTONES.md#{r['id'].lower()}) | {r['title']} | {earned(r)} / {r['points']} | {r['status']} |" for r in data['milestones']]
    card += ['', '## Execution queue', '',
             'An unmet dependency blocks acceptance, not exploratory work. The queue is a default ordering; a new theorem or the user’s direction can change it. Check the concrete prerequisite before starting.', '',
             '| Priority | Milestone | Unmet dependencies |', '|---:|---|---|']
    for i, identifier in enumerate(data['priority_queue'],1):
        row = by_id[identifier]
        blocked = [d for d in row['dependencies'] if by_id[d]['status'] != 'met']
        card.append(f"| {i} | [{identifier} · {row['title']}](MILESTONES.md#{identifier.lower()}) | {', '.join(blocked) or 'None'} |")
    card += ['', '## Checks actually run during this assessment', '']
    for check in data['fresh_checks_this_assessment']:
        card += [f"- `{check['command']}` — {check['result']} Limit: {check['limit']}"]
    card += ['', 'The score validator additionally checks totals, unique IDs, evidence references, dependency cycles, met prerequisites and generated-file consistency. It does not replace the mathematical acceptance tests.', '']

    detail = prelude + ['# Sixty research milestones', '', stamp, '',
        'Each gate is binary: every listed acceptance condition must hold to earn its full weight. “Partial” records useful existing work but earns zero. Evidence links are pinned; source descriptions state their exact limits.', '',
        'See [scope](SCOPE.md), [current score](SCORECARD.md), [structural assessment](CURRENT_STAGE.md) and [execution protocol](AGENTS.md).', '']
    for track in data['tracks']:
        detail += [f"## {track['id']}. {track['title']} — {track['maximum']} points", '']
        for row in data['milestones']:
            if row['track'] != track['id']:
                continue
            detail += [f"### {row['id']}", '', f"**{row['title']}** — {earned(row)}/{row['points']} points; **{row['status']}**.", '',
                f"**Scope:** {row['scope']}", '', '**Completion test:**', '']
            detail += [f'{i}. {s}' for i,s in enumerate(row['acceptance'],1)]
            detail += ['', f"**Current assessment:** {row['current_assessment']}", '',
                       '**Evidence:** ' + '; '.join(f"[{e}]({evidence[e]['url']})" for e in row['evidence']) + '.', '',
                       '**Prerequisites:** ' + (', '.join(f'[{d}](#{d.lower()})' for d in row['dependencies']) or 'None') + '.', '',
                       f"**Next calculation:** {row['next_action']}", '',
                       f"**Required artifact:** {row['required_artifact']}", '']

    sources = prelude + ['# Evidence register', '', stamp, '',
        'This register covers the scored assessment. A07 requires a larger claim ledger covering every advertised result in the repository and its linked public commits. Inclusion here alone does not certify a claim.', '']
    for key, item in evidence.items():
        sources += [f'## {key}', '', f"[{item['title']}]({item['url']})", '', item['finding'], '', f"**Review scope:** {item['review']}", '']
    return dict(zip(GENERATED, ('\n'.join(card), '\n'.join(detail), '\n'.join(sources))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true', help='regenerate the three rendered Markdown files')
    parser.add_argument('--ready', action='store_true', help='also list unearned gates whose declared prerequisites are met')
    args = parser.parse_args()
    data = json.loads((HERE/'milestones.json').read_text())
    evidence = json.loads((HERE/'evidence.json').read_text())
    scope = json.loads((HERE/'scope.json').read_text())
    if scope['scope_version'] != data['rubric_version']:
        raise ValueError('scope and rubric versions disagree')
    groups = scope['core_full_jet_groups']
    if not groups or len(groups) != len(set(groups)):
        raise ValueError('empty or duplicated gauge-group coverage')
    if scope['projection_constructions_allowed'] or scope['new_su2_sp1_so3_scans_allowed']:
        raise ValueError('scope conflicts with the current research constraints')
    if scope['coverage']['empty_denominator_is_success']:
        raise ValueError('empty coverage denominators cannot count as success')
    by_id = validate(data, evidence)
    rendered = render(data, evidence)
    for filename, content in rendered.items():
        path = HERE/filename
        if args.write:
            path.write_text(content)
        elif not path.exists() or path.read_text() != content:
            raise ValueError(f'{filename} is stale; run python3 research-roadmap/score.py --write')
    summary = totals(data)
    summary.update(rubric_version=data['rubric_version'], baseline_commit=data['baseline_commit'],
                   validation='passed', generated_files='written' if args.write else 'matched')
    if args.ready:
        summary['ready'] = [r['id'] for r in data['milestones'] if r['status'] != 'met' and
                            all(by_id[d]['status']=='met' for d in r['dependencies'])]
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (KeyError, TypeError, ValueError) as error:
        print(f'Roadmap validation failed: {error}', file=sys.stderr)
        sys.exit(1)
