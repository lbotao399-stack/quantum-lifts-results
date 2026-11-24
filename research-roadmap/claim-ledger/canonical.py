#!/usr/bin/env python3
"""Curated semantic records for the public undeformed-N4 claim ledger.

This file contains decisions, not a claim-mining heuristic.  Exact quotations
and unreviewed source spans live in the separate advertisement audit.  A
missing field remains missing; a link alone never upgrades a proof status.
"""
from fractions import Fraction

BASELINE = '19cd8718b752b04104d6ff0becee42d580b27ca3'
REPOSITORY = 'lbotao399-stack/quantum-lifts-results'


def reference(path, commit=BASELINE, anchor=None):
    result = {'kind': 'git', 'commit': commit, 'path': path,
              'url': 'https://github.com/' + REPOSITORY + '/blob/' + commit + '/' + path}
    if anchor:
        result['anchor'] = anchor
        result['url'] += '#' + anchor
    return result


def native(theta, dotted, length, **extra):
    total = sum(theta)
    return {'coordinate_order': ['theta1', 'theta2', 'theta3', 'n_plus', 'n_minus'],
            'native_weight': [*theta, *dotted], 'letter_length': length,
            'twice_energy': total + length + 2 * sum(dotted),
            'twice_R': [2*t - total + length for t in theta],
            'twice_J_left': total - length + sum(dotted),
            'twice_J_right': dotted[0] - dotted[1],
            'physical_level': 2 * total + 3 * sum(dotted), **extra}


def power(number):
    return native([0, 0, 0], [number, number], 2 * number)


def incoming(target, order):
    t1, t2, t3, p, m = target['native_weight']
    return native([t1 + order, t2 + order, t3 + order],
                  [p - order, m - order], target['letter_length'] + order - 1)


def outgoing(source, order):
    t1, t2, t3, p, m = source['native_weight']
    return native([t1 - order, t2 - order, t3 - order],
                  [p + order, m + order], source['letter_length'] + 1 - order)


def explicit_scope(text, completeness='stated complete block'):
    return {'description': text, 'completeness': completeness,
            'physical_identification': 'not inferred from a formal result'}


def claims():
    result = []

    def add(identifier, group, title, proposition, value, charges, order,
            certificate, *, direction='not_applicable', source=None, target=None,
            state='current', proof='published_certificate', field='Q',
            normalization='Dimension/rank claim: unchanged by nonzero basis rescaling.',
            scope=None, supersedes=(), limitations=(), unresolved=(), advertisements=()):
        certificates = certificate if isinstance(certificate, list) else [certificate]
        result.append({'id': identifier, 'theory': 'undeformed four-dimensional N=4 SYM',
                       'group': group, 'title': title,
                       'proposition': {'key': proposition, 'value': value},
                       'charges': charges,
                       'source': source if source is not None else {'kind': 'not_applicable', 'reason': 'This record is not a map claim.'},
                       'target': target if target is not None else {'kind': 'not_applicable', 'reason': 'This record is not a map claim.'},
                       'order': order, 'direction': direction, 'field': field,
                       'normalization': normalization,
                       'scope': scope or explicit_scope(title),
                       'certificate': certificates, 'supersession_state': state,
                       'proof_status': proof, 'supersedes': list(supersedes),
                       'limitations': list(limitations), 'unresolved_metadata': list(unresolved),
                       'advertisements': list(advertisements)})

    supplement = reference('supplement/supplement.tex')
    discussion = reference('RESEARCH_DISCUSSION.md')
    # Historical frontier is a finite-field campaign claim.  The absent
    # per-sector receipt is not silently replaced by a complete modern census.
    for group, bound, partial in [
            ('SU(2)', 33, '35 by restriction detector; chains at 39 and 40'),
            ('SU(3)', 30, 'detector through 34'), ('SU(4)', 28, 'none stated'),
            ('SU(5)', 28, '29:26/29; 30:3/53'), ('SU(6)', 30, '31:1/38'),
            ('SU(7)', 30, '31:19/37; 32:14/62'),
            ('SU(8)', 31, '32:45/56; 33:19/45'),
            ('SU(9)', 33, '34:41/64; 35:18/50'),
            ('SU(10)', 34, '35:33/38; 36:39/72'),
            *[('SU(%d)' % n, 36, 'none stated') for n in range(11, 16)]]:
        add('historical-frontier-' + group.lower().replace('(', '').replace(')', ''),
            group, 'Historical complete frontier through level %d' % bound,
            'historical-census/' + group, {'inclusive_level': bound, 'partial_levels': partial},
            {'domain': 'Every sector advertised in the historical census through physical level %d' % bound},
            'classical', supplement, state='historical', proof='finite_field_campaign',
            field='F_2147483647; exceptions retain their own stated fields',
            scope=explicit_scope('Supplement, sec-E-SU-frontier; this is a dated claimed frontier, not the new level-36 coverage denominator.'),
            unresolved=['The table does not identify an immutable complete per-sector input/receipt for every advertised sector.'],
            advertisements=['supplement:sec-E-SU-frontier', 'README:Theories'])

    # The explicit charge table is retained in its original normalization.
    classical_rows = [
        ('su2-o0','SU(2)','O_0',[4,4,4],[0,0],7,{'class_dimension':1},'Q'),
        ('su2-o1','SU(2)','O_1',[6,6,6],[0,0],9,{'named_class':'O_1'},'F_2147483647,F_2147483629'),
        ('su2-p39','SU(2)','P39',[5,5,5],[3,0],9,{'cohomology':13,'known_products_descendants':12,'remaining':1},'F_2147483647,F_2147483629'),
        ('su2-x40','SU(2)','X40',[7,5,5],[1,1],10,{'cohomology':8,'Cartan_rank':8,'known_products':7},'F_2147483647,F_2147483629'),
        ('su3-x7','SU(3)','X7',[2,2,2],[2,2],7,{'cohomology':26,'Cartan_rank':26,'known_products':25,'highest_weight_blocks':35},'F_2147483647'),
        ('su3-x7-extreme','SU(3)','X7 pure-jet member',[0,0,0],[5,5],7,{'trace_generators':26,'nonzero_class':True},'Q'),
        ('su3-primary-54411','SU(3)','Primary (5,4,4;1,1) length 10',[5,4,4],[1,1],10,{'cohomology':2,'Cartan_rank':0},'F_2147483647'),
        ('su3-primary-74410','SU(3)','Primary (7,4,4;1,0) length 11',[7,4,4],[1,0],11,{'cohomology':2,'beyond_F1_descendant':1,'Cartan_rank':0},'F_2147483647'),
        ('su3-primary-64411','SU(3)','Primary (6,4,4;1,1) length 11',[6,4,4],[1,1],11,{'cohomology':3,'Cartan_rank':0},'F_2147483647'),
        ('su3-primary-75410','SU(3)','Primary (7,5,4;1,0) length 11',[7,5,4],[1,0],11,{'cohomology':4,'Cartan_rank':0},'F_2147483647'),
        ('su4-threshold','SU(4)','First full-jet SU4 class',[4,2,2],[2,2],9,{'fortuitous_dimension':1},'Q'),
        ('sp2-threshold','Sp(2)=USp(4)','Sp2 threshold',[3,2,2],[1,1],7,{'fortuitous_dimension':1,'Cartan_rank':0},'F_2147483647'),
        ('so7-of','SO(7)','O_f',[3,3,3],[0,0],8,{'cohomology':1,'Cartan_rank':1},'F_2147483647'),
        ('sp3-threshold','Sp(3)=USp(6)','Sp3 threshold',[5,2,2],[1,1],9,{'cohomology':1,'fortuitous_dimension':1},'F_2147483647'),
        ('sp3-omega','Sp(3)=USp(6)','Known omega',[5,0,0],[3,3],9,{'cohomology':1,'Cartan_rank':0},'F_2147483647'),
        ('so9-threshold','SO(9)','SO9 threshold',[5,3,3],[0,0],10,{'cohomology':6,'beyond_Sp4':1},'F_2147483647'),
        ('so11-threshold','SO(11)','SO11 threshold',[7,3,3],[0,0],12,{'cohomology':23,'fortuitous_dimension':1},'F_2147483647')]
    for identifier, group, name, theta, dotted, length, value, field in classical_rows:
        add('historical-class-' + identifier, group, name + ' classical claim',
            'classical/' + identifier, value, native(theta,dotted,length),
            'classical', supplement, state='historical', field=field,
            proof='published_statement' if field == 'Q' else 'finite_field_campaign',
            unresolved=['Historical table needs its complete native input and certificate receipt linked at the row level.'],
            advertisements=['supplement:sec-E-SU-classes' if group.startswith('SU') else 'supplement:sec-E-SOSp-first', 'README:Theories'])
    # Exact later certificates supersede only the claim they actually decide.
    omega = native([5,0,0],[3,3],9,
                   original_package_order=['n_plus','n_minus','theta1','theta2','theta3'],
                   original_package_tuple=[3,3,5,0,0])
    omega_ref=reference('algorithms/fortuitous_duals_20261123/README.md')
    add('sp3-omega-exact-fortuity','Sp(3)=USp(6)',
        'Complete known-omega finite-colour relation quotient','classical/sp3-omega',
        {'cohomology':1,'fortuitous_dimension':1,'stable_source':2870,
         'finite_source':2689,'source_relations':181,'relation_image_rank':116,
         'stable_outgoing_rank':1198,'finite_outgoing_rank':1081,'stable_cohomology':0,
         'detector_pairing':'1'},omega,'classical',omega_ref,
        normalization='Original 4264-term omega; the published forty-entry detector pairs to one.',
        supersedes=['historical-class-sp3-omega'],
        limitations=['No new class or generator is created by certifying the already known omega.'],
        advertisements=['discussion:DUAL-CHATGPT-20261123'])
    add('sp-rank-window','Sp(N), N>=1','Known-omega grading vanishes for every rank at least four',
        'classical/sp-omega-rank-family',{'dimensions_ranks_1_2_3_4':[0,0,1,0],'all_ranks_at_least_4':0},
        omega,'classical',reference('algorithms/fortuitous_duals_20261123/rank_persistence/README.md'),
        scope=explicit_scope('One fixed complete native grading. Source injectivity at rank four and symplectic block embedding prove the unbounded-rank exclusion.'),
        advertisements=['discussion:DUAL-CHATGPT-20261123'])

    for identifier, group, theta, dotted, length, order, rank in [
        ('so7-of','SO(7)',[3,3,3],[0,0],8,1,1),
        ('so7-433','SO(7)',[4,3,3],[0,0],8,1,1),
        ('so7-mixed311','SO(7)',[3,1,1],[2,2],8,1,1),
        ('so9-threshold','SO(9)',[5,3,3],[0,0],10,1,1),
        ('so11-threshold','SO(11)',[7,3,3],[0,0],12,1,1)]:
        source=native(theta,dotted,length)
        add('historical-lift-'+identifier,group,'Historical first-order '+identifier,
            'map/'+identifier+'/outgoing/1',{'rank':rank},source,order,supplement,
            direction='outgoing',source=source,target=outgoing(source,order),state='historical',
            field='Q' if identifier=='so7-mixed311' else 'F_2147483647; SO7 O_f also checked at two other primes',
            proof='published_statement' if identifier=='so7-mixed311' else 'finite_field_campaign',
            unresolved=['The historical table lacks the exact row-specific input/certificate pin.'],
            advertisements=['supplement:sec-E-SOSp-lift','README:Theories'])
    for identifier,group,theta,dotted,length,field in [
        ('su2-p39','SU(2)',[5,5,5],[3,0],9,'F_2147483647,F_2147483629'),
        ('su2-x40','SU(2)',[7,5,5],[1,1],10,'F_2147483647; two primes only through order two'),
        ('su3-x7','SU(3)',[2,2,2],[2,2],7,'F_2147483647'),
        ('su4-threshold','SU(4)',[4,2,2],[2,2],9,'F_2147483647'),
        ('sp2-threshold','Sp(2)=USp(4)',[3,2,2],[1,1],7,'F_2147483647'),
        ('sp3-threshold','Sp(3)=USp(6)',[5,2,2],[1,1],9,'F_2147483647')]:
        charge=native(theta,dotted,length)
        add('historical-survival-'+identifier,group,'Historical advertised all-order survival of '+identifier,
            'survival/'+identifier,{'all_perturbative_orders':True},charge,'all_perturbative_orders',
            supplement,source=charge,target=charge,state='historical',field=field,
            proof='finite_field_campaign',
            scope=explicit_scope('The historical table advertises a formal all-order result; this record retains its finite-field certification level.'),
            unresolved=['No complete rational replay was located for this historical table entry.'],
            advertisements=['supplement:sec-E-SU-lift' if group.startswith('SU') else 'supplement:sec-E-SOSp-first','README:Theories'])
    add('historical-su3-two-class-partial-survival','SU(3)',
        'Two classes in (5,4,4;1,1), length ten, survive through order three',
        'survival/su3-primary-54411/through3',{'incoming_rank':0,'outgoing_rank':0},
        native([5,4,4],[1,1],10),[1,2,3],supplement,
        direction='both',source=native([5,4,4],[1,1],10),target=native([5,4,4],[1,1],10),
        field='F_2147483647',state='historical',proof='finite_field_campaign',
        limitations=['README labels the fourth-order continuation conditional; it is not an unconditional rational theorem.'],
        unresolved=['The complete historical rational chain and row-specific certificate are not attached.'])

    # Current large certificates, retaining formal versus physical scope.
    r007=native([7,4,4],[0,0],10,charge_conjugation='odd')
    add('su3-r007-formal-survival','SU(3)','R-007 has a non-torsion rational formal completion',
        'survival/su3-r007',{'all_perturbative_orders':True,'formal_torsion':False},
        r007,'all_perturbative_orders',reference('nonperturbative-lift/contributions/dot-r007-all-order/ALL_ORDER_PROOF.md'),
        direction='both',source=r007,target=r007,field='Q[[formal perturbative parameter]]',
        normalization='The frozen 6139-term R-007 representative; block-dual pairing 31104.',
        limitations=['Existence does not print all repaired coefficients and does not determine a physical instanton matrix.'],
        advertisements=['discussion:N4-NP-DOT-006','discussion:N4-NP-CODEX-004'])
    f1={'parameter':'integer n>=3','native_weight':['n+4',4,4,0,0],
        'letter_length':'n+7','physical_level':'24+2*n',
        'twice_energy':'2*n+19','twice_R':['2*n+3',3,3],
        'twice_J_left':5,'twice_J_right':0,
        'seed_records':[native([7,4,4],[0,0],10),native([8,4,4],[0,0],11)]}
    add('su3-f1-formal-tower','SU(3)','F1 has a non-torsion formal survivor at every stated level',
        'survival/su3-f1-family',{'all_n_at_least_3':True},f1,'all_perturbative_orders',
        reference('nonperturbative-lift/contributions/dot-su3-f1-tower/ALL_ORDER_PROOF.md'),
        direction='both',source=f1,target=f1,field='Q[[formal perturbative parameter]]',
        normalization='The separately pinned 6139-term level-30 and 8640-term level-32 seeds and their marked scalar products.',
        limitations=['This is an existence and non-torsion theorem, not physical all-coupling protection.'],
        advertisements=['discussion:N4-NP-DOT-007','README:Theories'])
    b1={'parameter':'integer n>=3','native_weight':['n+5',5,4,0,0],
        'letter_length':'n+8','physical_level':'28+2*n',
        'twice_energy':'2*n+22','twice_R':['2*n+4',4,2],
        'twice_J_left':6,'twice_J_right':0,
        'seed_records':[native([8,5,4],[0,0],11),native([9,5,4],[0,0],12)]}
    add('su3-b1-formal-tower','SU(3)','B1 has a non-torsion formal survivor at every stated level',
        'survival/su3-b1-family',{'all_n_at_least_3':True,'multiplier_rank':53},b1,
        'all_perturbative_orders',reference('algorithms/su3_j34_20261118/local_relations/README.md','c7e696aca2745db5283a0be1308b6f73d73797be'),
        direction='both',source=b1,target=b1,field='Q[[formal perturbative parameter]]',
        normalization='Pinned level-34 and level-36 seeds; pairings 120294720 and 19780588800; scalar dressing multiplier 6^k.',
        limitations=['Physical primary matching and positive dotted-grade instanton partners remain open.'],
        advertisements=['discussion:N4-NP-CODEX-004','README:Theories'])
    nine=native([5,2,2],[2,2],9)
    add('su3-nine-letter-formal-survival','SU(3)','The exceptional nine-letter target survives every perturbative page',
        'survival/su3-nine-letter',{'all_perturbative_orders':True},nine,
        'all_perturbative_orders',reference('algorithms/su3_j30_classical_ranks/bmn_return/README.md','acf631f4f92c2c271bbf1d6076d5d576e5dd0b7b'),
        direction='both',source=nine,target=nine,field='Q[[formal perturbative parameter]]',
        normalization='34492-term representative; four reverse Taylor operations give 240 times a nonzero 477-term BMN class.',
        advertisements=['discussion:N4-NP-CODEX-REVERSE'])
    add('su3-complete-second-source-zero','SU(3)','Complete source block gives zero incoming second differential to the nine-letter target',
        'map/su3-nine-letter/incoming/2',{'rank':0,'even_source_cohomology':4,'odd_source_cohomology':1},
        nine,2,[reference('algorithms/su3_j30_classical_ranks/bmn_return/README.md','c7e696aca2745db5283a0be1308b6f73d73797be'),
                reference('algorithms/su3_j30_classical_ranks/bmn_return/source-products.json','c7e696aca2745db5283a0be1308b6f73d73797be'),
                reference('algorithms/su3_j30_classical_ranks/bmn_return/verify.py','c7e696aca2745db5283a0be1308b6f73d73797be')],
        direction='incoming',source=incoming(nine,2),target=nine,
        limitations=['All five source classes are represented by four graviton products and the completed odd R-007.',
                     'README wording "At lengths nine and ten, respectively" describes adjacent map ranks ambiguously. check_source in the immutable verifier and the product input both explicitly fix the cohomology source length to ten in both colour parities.'],
        advertisements=['discussion:N4-NP-CODEX-004'])
    seven=native([3,0,0],[4,4],7,charge_conjugation='odd',dotted_spin=0)
    add('su3-seven-letter-generator-count','SU(3)','One length-seven generator remains after all proper products',
        'classical/su3-seven-letter-generator',
        {'cohomology':34,'proper_product_dimension':33,'generator_quotient_dimension':1,
         'cochain_dimensions':[1169,1896,1588],'incoming_rank':825,'outgoing_rank':1037,
         'factor_pairs':56,'factor_sectors':112,'factor_matrices':114},seven,'classical',
        reference('algorithms/su3_j30_classical_ranks/instanton_probe/outgoing_channels.md','a90728a6ca04d74f5475adae4bf37c8d94f34ac0'),
        limitations=['The count alone did not print the remaining generator; the 2026-11-24 packet supplies that separate object.'],
        advertisements=['discussion:SU(3) level 30, one full-jet ring generator'])
    add('su3-seven-letter-generator-explicit','SU(3)','Exact representative and full dual for the remaining seven-letter generator',
        'explicit/su3-seven-letter-generator',
        {'representative_terms':12224,'dual_entries':349,'dual_pairing':'1',
         'fortuitous_quotient':'2311-2310=1','proper_product_dimension':33},seven,'classical',
        [reference('algorithms/su3_generator_20261124/PROOF.md','eac7db42a6e27f4900826fe7b42cc75365698dc3'),
         reference('algorithms/su3_generator_20261124/MANIFEST.json','eac7db42a6e27f4900826fe7b42cc75365698dc3'),
         reference('algorithms/su3_generator_20261124/reproduce.py','eac7db42a6e27f4900826fe7b42cc75365698dc3'),
         reference('algorithms/su3_generator_20261124/independent_verify.py','eac7db42a6e27f4900826fe7b42cc75365698dc3')],
        proof='accepted_exact_certificate',
        normalization='Full rational representative normalized by the complete rational dual to pairing one.',
        limitations=['Physical superconformal parent and instanton role are separate questions.'])

    add('g2-level30-column','G2','The complete level-30 (5,5,5) column is exhausted by known products',
        'classical/g2-j30-555',{'dimensions_by_spin_0_through_9':[0,0,0,0,0,4,0,0,0,0],
                              'middle_dimension':3030,'outgoing_rank':1960,'incoming_rank':1066,
                              'known_product_augmented_rank':1070},
        {'native_weight':[5,5,5,0,0],'lengths_by_twice_spin':{'0':15,'1':14,'2':13,'3':12,'4':11,'5':10,'6':9,'7':8,'8':7,'9':6},
         'nonzero_block':native([5,5,5],[0,0],10),'physical_level':30},
        'classical',reference('algorithms/bmn_rank_certificates_20261117/README.md'),
        limitations=['Some required source-minor entries are retained production evaluations; C03 literal-source reconstruction is still blocked.'],
        advertisements=['README:November 18 BMN certificate release','STATUS:Verified BMN phase'])
    add('bmn-j30-index-screen',['Sp(3)=USp(6)','G2'],'Forty-one product sectors per group pass the level-30 index screen',
        'screen/bmn-j30',{'product_sectors_per_group':41,'weight_columns':12,'negative_bounds':0},
        {'domain':'The 41 sectors and 12 native weight columns in the pinned CSV.',
         'input':reference('data/N4_Sp3_G2_BMN_j30_direct_index_20261117.csv')},
        'classical',reference('algorithms/bmn_rank_certificates_20261117/README.md'),
        proof='finite_screen',scope=explicit_scope('Index minus certified product-rank bounds; cancelling classes remain unresolved.','complete stated screen, not complete cohomology'),
        advertisements=['STATUS:Verified BMN phase'])

    orth=reference('reviews/dot-orthogonal-d3/README.md','e1d184d388c588c0a5a825a19ecbca4e24427242')
    for group,num,value in [('SO(11)',6,'3072000'),('SO(12)',7,'-33177600'),('SO(13)',7,'-921600')]:
        charge=power(num)
        add(group.lower().replace('(','').replace(')','')+'-power%d-d3'%num,
            group,'Third differential into the defining-trace glueball power %d'%num,
            'map/'+group+'/power%d/incoming/3'%num,{'nonzero':True,'coefficient':value},charge,3,orth,
            direction='incoming',source=incoming(charge,3),target=charge,
            normalization='Dot normalization in the pinned full chain; defining representation trace.',
            scope=explicit_scope('Complete stated earlier-page certificate and full chain; SO12 source is reflection even.'),
            advertisements=['discussion:D003','discussion:Scope and current results'])
    add('so11-power6-d3-codex','SO(11)','Independent Codex source hits the same SO11 target',
        'map/SO(11)/power6/incoming/3/codex-normalization',{'nonzero':True,'coefficient':'8640'},power(6),3,
        reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/so11_j36/README.md','a90728a6ca04d74f5475adae4bf37c8d94f34ac0'),
        direction='incoming',source=incoming(power(6),3),target=power(6),
        normalization='Codex normalization: -17280+25920=8640. Coupling relation hbar_D=-2*hbar_C; not a second target.',
        limitations=['Different source normalizations and representatives do not multiply target or generator counts.'],
        advertisements=['discussion:D006','discussion:SU(3) level 30, one full-jet ring generator'])
    add('so13-power8-d3-product','SO(13)','The eighth-power third differential is the seventh-power product consequence',
        'map/SO(13)/power8/incoming/3',{'nonzero':True,'coefficient':'-921600','independent_generator_increment':0},
        power(8),3,reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/so13_j48/README.md','cc1346ec885eb45df78c31475fbe57155346fc38'),
        direction='incoming',source=incoming(power(8),3),target=power(8),
        normalization='Dot seventh-power chain multiplied by the same defining-trace glueball.',
        advertisements=['discussion:D008','discussion:Scope and current results'])
    add('so13-power8-earlier-pages','SO(13)','The eighth power survives the first and second incoming maps',
        'map/SO(13)/power8/incoming/1-2',{'rank':0,'complete_second_relations':301},power(8),[1,2],
        reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/so13_j48/README.md','cc1346ec885eb45df78c31475fbe57155346fc38'),
        direction='incoming',source=[incoming(power(8),1),incoming(power(8),2)],target=power(8),
        normalization='Second dual corrected by coefficients -9/1280 and -1/80 in the two specified distinct-plane coordinates.',
        advertisements=['discussion:D007','discussion:SO13-48-local','discussion:D008'])
    mixed_parent=native([7,5,5],[0,0],16)
    mixed_source=native([5,3,3],[3,3],16)
    mixed_target=native([2,0,0],[6,6],14)
    mixed_ref=reference('reviews/orthogonal-mixed-lifts/bmn13/THIRD.md','21ca40461b1bcacdf591f5f2fd85d85a5db6f8f0')
    add('so13-mixed-parent-classical-quotient','SO(13)','Two classical degree-17 relation classes leave one third-page source',
        'classical/SO13/mixed-parent',{'relation_dimensions_lengths_13_to_17':[0,9,35,69,42],
         'closed_relation_rank_length16':29,'relation_image_rank_length17':40,
         'relation_quotient_dimension':2,'first_differential_rank':1,'second_differential_rank':0},
        mixed_parent,[0,1,2],reference('reviews/orthogonal-mixed-lifts/bmn13/README.md','21ca40461b1bcacdf591f5f2fd85d85a5db6f8f0'),
        direction='outgoing',source=mixed_parent,target=[outgoing(mixed_parent,1),outgoing(mixed_parent,2)],
        normalization='Ordered two-source combination 358475625*O0(seed0)-832*O0(seed1). First evaluations 872415232 and 375888936960000.',
        scope=explicit_scope('Complete stated relation quotient modulo stable cycles; the retained combination is a finite source class.'),
        advertisements=['issue:7075528933'])
    add('so13-mixed-level40-d3','SO(13)','A separate mixed level-40 third differential has exact nonzero pairing',
        'map/SO13/mixed-level40/incoming/3',{'nonzero':True,'pairing':'-50677394964480000',
         'residues_1000003_1000033':[248799,367978],'candidate_parents':42298},mixed_target,3,mixed_ref,
        direction='incoming',source=mixed_source,target=mixed_target,
        normalization='Defining trace, native ordered-pair kernel K=1; target covector (6,-1,0); source is the six ordered upward substitutions of 358475625*O0(seed0)-832*O0(seed1).',
        limitations=['The same upward identity detects the j34 parent obstruction. This does not by itself transport Dot\'s separate j42 map downward.','No new replay of this public certificate is asserted by the ledger.'],
        advertisements=['issue:7076191991'])
    add('so13-mixed-parent-d3','SO(13)','The derivative-free level-34 parent has a nonzero third obstruction',
        'map/SO13/mixed-parent/outgoing/3',{'nonzero':True,'detected_after_six_upward_substitutions':'-50677394964480000'},mixed_parent,3,mixed_ref,
        direction='outgoing',source=mixed_parent,target=outgoing(mixed_parent,3),
        normalization='The same ordered source combination and verified upward substitution identities as the level-40 descendant.',
        limitations=['Parent and descendant do not constitute two independent generator families.'],
        advertisements=['issue:7076191991'])
    add('so13-mixed-kernel-exhaustion-conditional','SO(13)','Two mixed third images would exhaust the Cartan kernel under the downward Ward hypothesis',
        'conditional/SO13/mixed-target-kernel-exhaustion',{'terminal_dimension':3,'Cartan_rank':1,'Cartan_kernel_dimension':2,
         'independent_target_vector_from_Dot':[7,42,0],'annihilating_covector':[6,-1,0]},mixed_target,3,mixed_ref,
        direction='incoming',source=[mixed_source,incoming(power(7),3)],target=mixed_target,
        proof='conditional_formal_identification',
        normalization='The level-40 row (6,-1,0) annihilates (7,42,0).',
        limitations=['Transport of Dot\'s pure-power differential downward requires filtered superconformal equivariance. The independently nonzero level-40 pairing does not require this downward hypothesis.'],
        advertisements=['issue:7076191991'])
    add('so11-aligned-source-fortuity-zero','SO(11)','The aligned leading source combination has zero fortuitous obstruction',
        'fortuity/so11-aligned-combination',{'fortuitous_obstruction':0,'finite_class_nonzero':'unresolved'},
        incoming(power(6),3),'classical',reference('algorithms/fortuitous_duals_20261123/README.md'),
        normalization='9 Y_D(-2 hbar_C)+25600 Y_C(hbar_C); quotient 9*(1/6)+(25600/25200)*(-189/128)=0.',
        limitations=['The two old coordinate probes distinguish representatives but do not annihilate the full relation image.','The finite leading class and a complete stable quantum lift remain open.'],
        advertisements=['discussion:SO11-RELATIONS-CHATGPT-20261123'])
    for group,num,level,page in [('SO(19)',10,60,5),('SO(25)',13,78,6),('SO(26)',14,84,6)]:
        add(group.lower().replace('(','').replace(')','')+'-finite-page-report',group,
            'Reported finite-page survival of glueball power %d'%num,
            'survival/'+group+'/power%d/page%d'%(num,page),{'survives_to_page':page},power(num),
            list(range(1,page)),discussion,direction='incoming',source=[incoming(power(num),r) for r in range(1,page)],target=power(num),
            proof='reported_unintegrated',unresolved=['A complete separately pinned native proof/replay packet for this finite-page statement has not been integrated.'],
            limitations=['Later allowed orders remain open.'],advertisements=['discussion:D014'])
    add('so14-power8-d3-report','SO(14)','Reported SO14 eighth-power third-order hit',
        'map/SO(14)/power8/incoming/3',{'nonzero':True},power(8),3,discussion,
        direction='incoming',source=incoming(power(8),3),target=power(8),
        proof='reported_unintegrated',normalization='D014 states that the signed odd-cycle contribution is included; the literal normalized witness is not attached here.',
        unresolved=['The exact public source, chain, coefficient normalization and full earlier-page certificate require a native immutable proof pin.'],
        advertisements=['discussion:D014'])
    add('orthogonal-first-invisible-conjecture','SO(N), N>=7','First Cartan-invisible pure glueball power is conjectured to lift at some finite order',
        'conjecture/orthogonal-first-invisible',{'conjectured':True},
        {'parameter':'integer N>=7','power':'floor(N/2)+1','native_weight':[0,0,0,'floor(N/2)+1','floor(N/2)+1'],
         'letter_length':'2*floor(N/2)+2','hypothesis':'The stated power is classically nonzero.'},
        'some finite positive perturbative order',discussion,proof='conjecture',
        scope=explicit_scope('D014 is an SO-only conjecture; finite samples are evidence, not a uniform theorem.','conjecture'),
        advertisements=['discussion:D014'])
    add('orthogonal-first-family','SO(N)','Uniform first-differential pure-glueball family',
        'map/orthogonal-first-family',
        {'coefficient':'3*2^(p-5)*(2*N-3*p-4)/(p-4)!',
         'source_closed_for':'p>=4 and N<p+5',
         'target_nonzero_for':'p<N-2',
         'all_N_at_least_9_first_hit_powers':['N-4','N-3'],
         'one_source_zero_endpoint':{'N':8,'p':4}},
        {'parameters':'integers p>=4 and N>=7; the source is closed for N<p+5; the target is nonzero for p<N-2',
         'source_native_weight':[1,1,1,'p-1','p-1'],'source_letter_length':'2*p',
         'target_native_weight':[0,0,0,'p','p'],'target_letter_length':'2*p',
         'target_twice_energy':'6*p','target_twice_R':['2*p','2*p','2*p'],
         'target_twice_J_left':0,'target_twice_J_right':0,'physical_level':'6*p',
         'source_twice_energy':'6*p-1','source_twice_R':['2*p-1','2*p-1','2*p-1'],
         'source_twice_J_left':1,'source_twice_J_right':0,
         'native_input_reference':reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/first_family.md','d320363f2eb54898e6ff554943b0d6d66705c120')},
        1,reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/first_family.md','d320363f2eb54898e6ff554943b0d6d66705c120'),
        direction='incoming',source={'native_weight':[1,1,1,'p-1','p-1'],'letter_length':'2*p',
                                     'recipe':'One scalar-marker determinant, six single ghost paths, and p-4 double ghost paths; the derivative is a degree-(p+5) polarized Pfaffian coefficient.'},
        target={'native_weight':[0,0,0,'p','p'],'letter_length':'2*p','representative':'(Tr_N(lambda_plus lambda_minus))^p'},
        normalization='Published odd-trace reversal and marker order +H,+E,+F,-H,-E,-F; coefficient 3*2^(p-5)*(2*N-3*p-4)/(p-4)!.',
        limitations=['The zero of this one source at (N,p)=(8,4) is not a complete zero theorem.'],
        advertisements=['discussion:D007','discussion:D014'])

    # Exceptional chronology: a withdrawal is not a full-source zero.
    f4first=reference('algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md','03fc9a5d8a68a1b70072242a0770335df1d862d3')
    add('f4-power5-formal-survival','F4','Fifth glueball power survives all perturbative pages',
        'survival/F4/power5',{'all_perturbative_orders':True,'second_source_cohomology':'141-80-61=0'},power(5),
        'all_perturbative_orders',f4first,direction='both',source=power(5),target=power(5),
        normalization='Native F4 colour normalization and defining-trace glueball in the pinned proof.',
        limitations=['Its physical instanton current matrix element remains open.'],advertisements=['discussion:F4-CODEX-FIFTH'])
    add('f4-power6-earlier-page','F4','Sixth power has a nonzero first-page Ward covector',
        'map/F4/power6/incoming/1',{'rank':0,'target_nonzero':True,'independent_SO9_pairing':'720'},power(6),1,
        reference('algorithms/dual_singlets_20261118/f4/j36_source/first_page/README.md','af0578d98b9182fcd92f288f862e033d47765293'),
        direction='incoming',source=incoming(power(6),1),target=power(6),
        normalization='Native dual pairings -15313094596800 and 720; all 763 differentiated-scalar and 6920 mixed second-ghost coordinates tested.',
        advertisements=['discussion:F4-CODEX-E2'])
    withdrawal=reference('reviews/f4-j36-eight-albert-d2/README.md','7fcfb7175703dd1ee6dcf4a71f2885c0d7dc8c73')
    add('f4-power6-original-witness-withdrawn','F4','The original sixth-power witness is classically exact and withdrawn',
        'witness/F4/power6/original/d2',{'source_class':0,'corrected_second_pairing':0,
         'primitive_terms':6417,'collected_source_words':10401},power(6),2,
        reference('algorithms/dual_singlets_20261118/f4/j36_source/README.md','ee9e06494a926fb88baeb6fd875cee36762f63a7'),
        direction='incoming',source=incoming(power(6),2),target=power(6),state='withdrawn',proof='withdrawal',
        normalization='Original archive conventions; restored branch metrics -4,+2 give 283511008535040000/7-283511008535040000/7=0.',
        limitations=['This withdrawal predates the distinct eight-Albert witness and must not be merged with it.','The attached seventh-power product consequence was withdrawn with this source.'],
        advertisements=['discussion:D001','discussion:D002','issue:7074126654','issue:7074234264'])
    add('f4-power7-original-product-withdrawn','F4','The seventh-power consequence of the original exact source is withdrawn',
        'witness/F4/power7/original-product/d2',{'nonzero_witness_valid':False},power(7),2,
        reference('algorithms/dual_singlets_20261118/f4/j36_source/README.md','ee9e06494a926fb88baeb6fd875cee36762f63a7'),
        direction='incoming',source=incoming(power(7),2),target=power(7),state='withdrawn',proof='withdrawal',
        normalization='The original sixth-power source multiplied by the same closed glueball.',
        limitations=['No full-source seventh-power zero follows from withdrawal of this product.'],
        advertisements=['discussion:D002','issue:7074126654'])
    add('f4-power6-eight-albert-withdrawn','F4','The claimed eight-Albert second-order witness is withdrawn',
        'witness/F4/power6/eight-albert/d2',{'nonzero_witness_valid':False},power(6),2,withdrawal,
        direction='incoming',source=incoming(power(6),2),target=power(6),state='withdrawn',
        proof='withdrawal',normalization='Canonical exterior incidence order, separately for both dotted species.',
        limitations=['Withdrawal of this selected witness alone does not decide the complete source.'],
        advertisements=['discussion:F₄ glueball tower and the eight-Albert sign','issue:7105180683','issue:7105901721'])
    add('f4-power6-eight-albert-corrected-modular','F4','Corrected eight-Albert augmented determinants vanish at both archived primes',
        'witness/F4/power6/eight-albert/corrected-minor',
        {'residues_1009_1013':[0,0],'raw_rational_score':'-299015516814300/49','normalized_residues':[44,92]},power(6),2,
        reference('algorithms/dual_singlets_20261118/f4/j36_source/eight_albert_review/README.md','53462d3102e77eea433207ab1d0fd7e09a8303f6'),
        direction='incoming',source={'kind':'selected_eight_Albert_sources','native_charge':incoming(power(6),2)},target=power(6),
        field='Q for raw graph; F_1009,F_1013 for corrected augmented minors',
        scope=explicit_scope('Only the specified crossed source and the fixed classical section.','selected source'),
        limitations=['Two zero modular residues alone do not prove a rational zero.'],
        advertisements=['discussion:F₄ glueball tower and the eight-Albert sign','issue:7105901721','issue:7106200345'])
    add('f4-power6-bound15','F4','Earlier upper bound fifteen for the selected complete Taylor source',
        'bound/F4/power6/second-source',{'upper_bound':15,'dimensions':[3720,703,3002]},
        incoming(power(6),2),2,reference('algorithms/dual_singlets_20261118/f4/j36_source/global_boundaries/README.md','c3688a7ad9d17f5c949d22ff96d098846c5698b4'),
        state='historical',proof='rational_rank_bound',
        limitations=['3720-703-3002=15 is an upper bound, not an exact cohomology count.'],
        advertisements=['discussion:SU(3) level 30, one full-jet ring generator'])
    add('f4-power6-outgoing-classical-rank','F4','The complete selected thirteen-letter source has outgoing classical rank 703',
        'classical/F4/power6/second-source-outgoing-rank',{'cochain_dimension':3720,'outgoing_rank':703,'kernel_dimension':3017,
         'minor_residues_1009_1013':[1005,734],'Cartan_upper_bound_witness':'1728*5184=8957952'},
        incoming(power(6),2),'classical',reference('reviews/independent-f4/sixth-glueball-classical-rank/README.md','3b30b8ff242977debd90581c12f4a767d7875c6c'),
        source=incoming(power(6),2),target=native([2,2,2],[4,4],14),direction='outgoing',
        normalization='542 trace rows, 158 four-cubic Albert rows and three six-cubic Albert rows in the pinned native order.',
        advertisements=['issue:7088927273'])
    add('f4-power6-classical-products-second-zero','F4','Every contributing positive-degree classical product has zero incoming second image',
        'map/F4/power6/classical-products/incoming/2',{'rank_on_product_span':0,'complement_boundary_rank':55,
         'minor_residues_1009_1013':[87,474]},incoming(power(6),2),2,
        reference('reviews/independent-f4/sixth-glueball-classical-products/README.md','78d170a26afc2a16cae43f123df86d5880be6196'),
        source={'native_charge':incoming(power(6),2),'subspace':'Complete span of positive-degree classical products in Taylor diagram (4,4,3,3).'},
        target=power(6),direction='incoming',
        scope=explicit_scope('Entire stated classical product span; indecomposable sources require the later complete-source proof.','product subspace'),
        advertisements=['issue:7092236097'])
    newest=reference('algorithms/dual_singlets_20261118/f4/j36_source/focused_powers/sixth_survival.md','c778cd7c3ceb72ad8082d791b2aee5f0df082164')
    add('f4-power6-full-second-zero-report','F4','Complete incoming second differential is zero over the rationals',
        'map/F4/power6/incoming/2',{'rank':0,'kernel_dimension':3017,'boundaries':3016,
         'zero_obstruction_product':1,'minor_residues_1009_1013':[449,60],
         'third_source_cohomology':'1862-1434-428=0'},power(6),2,newest,
        direction='incoming',source=incoming(power(6),2),target=power(6),
        proof='accepted_published_theorem_with_independent_closed_factor_audit',
        normalization='Native ordered 3017 minor; 833*737 mod 1009=449 and 615*346 mod 1013=60.',
        limitations=['The root accepted D03 after auditing the published complete-source theorem and a new direct closed-factor proof. The complete 3017 minor remains an explicitly attributed published prerequisite; no new permitted full-minor replay is claimed.','The separate all-page formal permanence assertion retains its own proof status.'],
        supersedes=['f4-power6-bound15'],advertisements=['issue:7112142083','issue:7113458598'])
    add('f4-power6-formal-survival-report','F4','Newest public sixth-power permanence theorem in the relative complex',
        'survival/F4/power6',{'all_perturbative_orders':True},power(6),'all_perturbative_orders',newest,
        direction='both',source=power(6),target=power(6),
        proof='published_pending_independent_gate_review',
        normalization='The same native sixth glueball and complete-source proof as the second-order zero.',
        limitations=['Formal perturbative permanence is not an instanton result.'],
        advertisements=['issue:7112142083','issue:7113458598'])
    tower=reference('algorithms/dual_singlets_20261118/f4/j36_source/eight_albert_review/README.md','53462d3102e77eea433207ab1d0fd7e09a8303f6')
    add('f4-power7-first-zero','F4','Seventh power survives the first differential',
        'map/F4/power7/incoming/1',{'rank':0,'dual_pairing':'-194893931232000'},power(7),1,tower,
        direction='incoming',source=incoming(power(7),1),target=power(7),
        normalization='68 rational Ward coefficients; 1117 differentiated-scalar and 18948 mixed-jet residuals vanish.',
        advertisements=['discussion:F₄ glueball tower and the eight-Albert sign','issue:7105901721'])
    add('f4-power8-first-hit','F4','Eighth power is hit by the first differential',
        'map/F4/power8/incoming/1',{'nonzero':True,'homogeneous_rank':91,'augmented_rank':92,
         'homogeneous_minor_1009_1013':[156,382],'augmented_minor_1009_1013':[375,29]},power(8),1,tower,
        direction='incoming',source=incoming(power(8),1),target=power(8),
        normalization='Native homogeneous and augmented ordered Ward minors; nonzero finite-field minor proves rational nonvanishing.',
        advertisements=['discussion:F₄ glueball tower and the eight-Albert sign','issue:7105901721'])
    add('f4-power7-selected-exchange-zero','F4','518 marked Albert recipes supply no second-source classes',
        'source/F4/power7/selected-triple-exchange',{'formal_recipes':518,'outgoing_rank':191,'closed_combinations':327,
         'mark_differences':259,'colour_identities':68},incoming(power(7),2),2,
        reference('algorithms/dual_singlets_20261118/f4/j36_source/focused_powers/triple_exchange/README.md','c778cd7c3ceb72ad8082d791b2aee5f0df082164'),
        source={'native_charge':incoming(power(7),2),'selection':'Up to three ordered Albert exchanges, both marks and 259 words.'},target=power(7),direction='incoming',
        scope=explicit_scope('Complete stated 518-recipe family; not the full seventh-power source.','selected family'),
        limitations=['Independent seventh-power incoming orders two, three and four remain open.'],
        advertisements=['issue:7113458598'])
    add('f4-power7-later-open','F4','Complete incoming maps into the seventh power remain open',
        'open/F4/power7/incoming/2-4',{'unresolved_orders':[2,3,4]},power(7),[2,3,4],
        reference('algorithms/dual_singlets_20261118/f4/j36_source/focused_powers/seventh_ward.md','c778cd7c3ceb72ad8082d791b2aee5f0df082164'),
        source=[incoming(power(7),r) for r in [2,3,4]],target=power(7),direction='incoming',proof='open',
        normalization='The same independently tested native first-page target.',
        advertisements=['issue:7113458598'])

    # Remaining named compulsory inputs, with native charges even when open.
    for group,num,level in [('SO(11)',5,35),('SO(13)',6,41)]:
        charge=native([1,0,0],[num+1,num],2*num+2)
        add(group.lower().replace('(','').replace(')','')+'-one-scalar-earlier-pages',group,
            'One-scalar target at level %d reaches the third page'%level,
            'survival/'+group+'/one-scalar/page3',{'survives_to_page':3},charge,[1,2],
            reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/one_scalar/README.md','7deb71acbdf828bd46edad5ddf503c60cbdfd58b'),
            direction='incoming',source=[incoming(charge,1),incoming(charge,2)],target=charge,
            normalization='Defining trace S^%d Tr(lambda_dot_plus phi1).'%num,
            limitations=['All later incoming orders remain separate open work.'],advertisements=['discussion:D009'])
    quad=native([2,0,0],[7,7],16)
    add('so13-level46-two-classes','SO(13)','Both quadratic level-46 classes reach the third page',
        'survival/SO(13)/quadratic-level46/page3',{'second_page_dimension':2,'incoming_second_rank':0,
         'commuting_plane_rank':20,'augmented_rank':20},quad,[1,2],
        reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/quadratic_so13/README.md','6a6f5a708ab1f73f00f8baafd9a133e824fe8238'),
        direction='incoming',source=[incoming(quad,1),incoming(quad,2)],target=quad,
        normalization='The two target covectors (2,1,0) and (18,0,1) in the native three trace coordinates.',
        limitations=['The known scalar product is hit at third order; the independent other class requires its remaining complete decision.'],
        advertisements=['discussion:D011','discussion:D012'])
    add('g2-power3-formal-report','G2','D014 reports all-page formal survival of the third glueball power',
        'survival/G2/power3',{'all_perturbative_orders':True},power(3),'all_perturbative_orders',discussion,
        direction='both',source=power(3),target=power(3),proof='reported_unintegrated',
        normalization='Defining seven-dimensional trace in the stated G2 pure-glueball quotient.',
        unresolved=['The distinct complete native formal proof packet underlying D014 has not been integrated.'],
        limitations=['The physical protected-product partner argument is a separate conditional statement.'],
        advertisements=['discussion:D014'])
    add('g2-current-products-conditional','G2','Two current products conditionally exhaust the physical incoming channel of S^3',
        'physical/G2/power3/incoming',{'classical_dimensions_upper_bounds':[0,0,2,0],'determinant':'-36864','conditional_instanton_rank':0},
        power(3),'instanton charge one and stated protected-current implication',
        reference('algorithms/G2_current_product_span_20261120.md','97aafa6b253c479ef6370844151af51134b2684c'),
        direction='incoming',source={'native_input_reference':reference('algorithms/G2_current_product_span_20261120.md','97aafa6b253c479ef6370844151af51134b2684c')},target=power(3),
        proof='conditional_physical',normalization='Cartan probe rows (-64,0),(16,576).',
        limitations=['Requires the specified exact physical stress-multiplet products with their classical leading terms; no instanton integral was evaluated.'],
        advertisements=['discussion:N4-NP-CODEX-009'])
    for length in [11,13]:
        charge=native([6,6,6],[0,0],length)
        add('g2-level36-length%d-open'%length,'G2','G2 level-36 length-%d complete source and maps remain open'%length,
            'open/G2/level36/length%d'%length,{'complete_source_cohomology':'unresolved','higher_maps':'unresolved'},charge,
            'all kinematically allowed orders',discussion,source=charge,target={'native_map_neighbours':'Generated in native-frontiers.json.'},
            direction='outgoing',proof='open',advertisements=['discussion:G2-CODEX-TAYLOR36','discussion:G2-CODEX-CROSS36'])
    add('g2-length15-second-target-zero','G2','Length-fifteen level-36 second target has zero classical multiplicity',
        'classical/G2/level36/length15/second-target',{'classical_target_multiplicity':0},native([6,6,6],[0,0],15),2,
        reference('algorithms/dual_singlets_20261118/bmn_primary/axis_blocks/determinant/cross_bracket.md','d4608bf0270b7692dff90f154241ac44db9538f0'),
        direction='outgoing',source=native([6,6,6],[0,0],15),target=outgoing(native([6,6,6],[0,0],15),2),
        proof='conditional_formal_identification',
        limitations=['The physical outgoing conclusion uses the stated filtered Taylor equivariance.','Lengths eleven and thirteen remain unresolved.'],
        advertisements=['discussion:G2-CODEX-CROSS36'])
    add('usp4-four-odd-generator-bound','Sp(2)=USp(4)','Index constraints require at least four odd generators by level 48',
        'bound/Sp(2)/odd-generators-through48',{'lower_bound':4,'indices_at44_48':[-3,-2]},
        {'physical_levels':[44,48],'native_input_reference':reference('algorithms/dual_singlets_20261118/bmn_primary/second_odd_generator.md','d4608bf0270b7692dff90f154241ac44db9538f0')},
        'classical',reference('algorithms/dual_singlets_20261118/bmn_primary/second_odd_generator.md','d4608bf0270b7692dff90f154241ac44db9538f0'),
        proof='index_lower_bound',limitations=['This does not construct four explicit new generators; explicit-family increment is zero.'],
        advertisements=['discussion:G2-CODEX-CROSS36'])
    add('g2-scalar-determinant','G2','The six-scalar Gram determinant is a non-Cartan rational class',
        'classical/G2/scalar-gram',{'cochains':4,'boundaries':3,'cohomology':1,'Cartan_rank':0},
        native([2,2,2],[0,0],6),'classical',supplement,
        normalization='det Tr_7(gamma^I gamma^J), equivalently the explicitly normalized epsilon contraction with factor 1/6.',
        advertisements=['supplement:sec-E-G2-class','README:Theories'])
    add('g2-scalar-determinant-survival','G2','The six-scalar Gram determinant survives all perturbative orders',
        'survival/G2/scalar-gram',{'all_perturbative_orders':True},native([2,2,2],[0,0],6),
        'all_perturbative_orders',supplement,direction='both',source=native([2,2,2],[0,0],6),target=native([2,2,2],[0,0],6),
        proof='published_statement',normalization='Same scalar determinant as the classical record.',
        unresolved=['The complete historical higher-map native certificate has not been linked separately from the supplement statement.'],
        advertisements=['supplement:sec-E-G2-class','README:Theories'])

    # Physical records keep a partial zero strictly smaller than a full Ward map.
    cartan_target=native([5,2,2],[2,2],9,charge_conjugation='odd',
                         representative_id='su3-cartangraviton-5080-target')
    add('su3-cartangraviton-5080-target','SU(3)','A literal 5080-term Cartan-visible classical target has a formal graviton completion',
        'classical/SU3/cartangraviton-5080-target',{'nonzero':True,'terms':5080,'rational_Cartan_terms':34104,
         'first_Cartan_coefficient':'432/5','literal_closure_and_symmetry_tests':14,'formal_non_torsion':True},cartan_target,
        'classical and all perturbative orders',reference('nonperturbative-lift/contributions/dot-r007-t2-primary/PROOF.md'),
        direction='both',source=cartan_target,target=cartan_target,field='Q and Q[[formal perturbative parameter]]',
        normalization='Frozen 47-product rational recipe, followed by the four ordered Taylor operations in the proof.',
        limitations=['This Cartan-visible graviton target and the 34492-term Cartan-null nine-letter target have the same charges but are different classes.','The classical Taylor lowering tests do not establish a renormalized physical primary at nonzero coupling.'],
        advertisements=['discussion:N4-NP-DOT-008'])
    probe=reference('algorithms/su3_j30_classical_ranks/instanton_probe/flavour_modes.md','6a6f5a708ab1f73f00f8baafd9a133e824fe8238')
    add('su3-r007-classical-instanton-profile-zero','SU(3)','Bare charge-one classical profiles vanish in the specified orientations',
        'partial-physical/SU3/R007/classical-profiles',{'separated_classical_profile':0,'full_instanton_matrix':'open'},
        {'source':r007,'target':cartan_target},'instanton charge one',probe,direction='outgoing',source=r007,target=cartan_target,
        proof='partial_physical_exclusion',normalization='Pinned 6139 source words and 5080 target words; the two R-flavour occupation tests are stated in the proof.',
        scope=explicit_scope('Bare zero-vev classical profiles and the specified endpoint-only topologies.','partial physical contribution'),
        limitations=['Fluctuations, current-connected graphs, contact terms, mixing and integration-boundary effects are not decided.'],
        advertisements=['discussion:N4-NP-CODEX-005','discussion:N4-NP-CODEX-006'])
    add('su3-r007-minimum-current-flux-zero','SU(3)','Minimum endpoint-only Poincare-current flux vanishes',
        'partial-physical/SU3/R007/minimum-current-flux',{'stated_topology_flux':0,'full_instanton_matrix':'open'},
        {'source':r007,'target':cartan_target},'instanton charge one',
        reference('algorithms/su3_j30_classical_ranks/instanton_probe/background_exchange.md','d320363f2eb54898e6ff554943b0d6d66705c120'),
        direction='outgoing',source=r007,target=cartan_target,proof='partial_physical_exclusion',
        normalization='Native endpoint circuits with both terms of the adjoint propagator in the cited equation.',
        scope=explicit_scope('Only the minimum endpoint-only topology at smooth instanton size.','partial physical contribution'),
        limitations=['No full matched microscopic instanton differential has been calculated.'],
        advertisements=['discussion:N4-NP-CODEX-007'])
    add('su3-r007-classical-collision-support','SU(3)','The reduced scalar support of R-007 is the pair-collision locus',
        'support/SU3/R007/classical-collision',{'covering_directions':6,'five_suffice':False,'bare_discriminant_killing_power':9},
        r007,'classical',reference('algorithms/su3_j30_classical_ranks/instanton_probe/collision_support.md','1028aa7d18c83e591ac513c0715ccfa032cc85a4'),
        normalization='The fixed bare source and complete block-dual polynomial 31104 i h_1^3.',
        limitations=['This is a support theorem for the specified class, not a physical instanton pairing.'],
        advertisements=['discussion:F4 global boundaries and SU(3) collision support'])
    add('su3-r007-completed-discriminant','SU(3)','Completed R-007 has the stated finite discriminant-clearing relation',
        'support/SU3/R007/formal-discriminant',{'finite_total_exponent_exists':True,'first_clearing_exponent_at_least':6,
         'remaining_formal_powers':[1,3]},r007,'all_perturbative_orders',
        reference('algorithms/su3_j30_classical_ranks/instanton_probe/completed_discriminant.md','f823c91793fd3b804784a4218038c3947e573086'),
        field='Q[[formal perturbative parameter]]',
        normalization='The pinned R-007 representative and first-scalar discriminant.',
        limitations=['The total clearing exponent is existential.','The physical consequence requires a matched instanton derivation annihilating the stated graviton classes.'],
        advertisements=['discussion:N4-NP-CODEX-R007-LOCALIZATION'])
    add('historical-n4-campaign-totals-20261118','all campaign groups','November 18 campaign totals are a historical snapshot',
        'campaign/2026-11-18',{'audited_d1_highest_weight_pairings':5,'audited_generator_families':3,
                              'certified_d2':0,'certified_d3':0,'certified_microscopic_BPS_instanton':0},
        {'domain':'The November 18 audited undeformed-N4 campaign, not all repository results.'},
        'campaign snapshot',reference('STATUS.md'),state='historical',proof='historical_count',
        scope=explicit_scope('Counts exclude later orthogonal third-order results and later first-image multiplets.'),
        unresolved=['Historical campaign identifiers for every counted object are not yet joined to the semantic ledger.'],
        advertisements=['README:November 18 BMN certificate release','STATUS:Verified BMN phase'])
    # A07 records the exact advertised statement and its evidence state. It
    # does not impose a new requirement to replay every historical rank. The
    # following records already have explicit native charges or a quantified
    # aggregate domain; unavailable proof payloads remain scientific limits.
    for item in result:
        identifier=item['id']
        if identifier.startswith('historical-'):
            item['proof_status']='attributed_historical_claim'
            item['limitations'].extend('Historical evidence limitation: '+x for x in item['unresolved_metadata'])
            item['unresolved_metadata']=[]
            item['semantic_resolution']='The pinned table is the advertised statement, with its exact native charge/aggregate domain and stated certification field. A fresh mathematical acceptance is not asserted.'
            if identifier.startswith('historical-frontier-'):
                bound=item['proposition']['value']['inclusive_level']
                item['charges'].update(native_domain={
                    'native_coordinate_order':['theta1','theta2','theta3','n_plus','n_minus'],
                    'degrees':'all nonnegative integer degrees',
                    'level_constraint':'2*(theta1+theta2+theta3)+3*(n_plus+n_minus)<=%d'%bound,
                    'length_constraint':'vacuum separately; max(theta1,theta2,theta3)<=letter_length<=theta1+theta2+theta3+n_plus+n_minus',
                    'symmetry_scope':'The historical census reports highest-weight blocks and their complete multiplets.'})
        if identifier in ['so14-power8-d3-report','so19-finite-page-report','so25-finite-page-report','so26-finite-page-report','g2-power3-formal-report']:
            item['limitations'].extend('Unintegrated advertised proof: '+x for x in item['unresolved_metadata'])
            item['unresolved_metadata']=[]
            item['proof_status']='attributed_public_claim_without_integrated_replay'
            item['normalization']='D014 defines the glueball by the defining-representation trace. The stated nonvanishing or vanishing is invariant under nonzero source and target rescaling. No literal numerical pairing coefficient is attributed to this report.'
            item['semantic_resolution']='The target and each complete allowed source charge are explicit. The source is existentially quantified in the reported map statement; an unprinted witness is a proof-availability limitation, not an invented representative.'
        if identifier=='g2-scalar-determinant-survival':
            item['limitations'].extend('Historical evidence limitation: '+x for x in item['unresolved_metadata'])
            item['unresolved_metadata']=[]
            item['proof_status']='attributed_historical_claim'
            item['field']='Q for the classical determinant; the higher-map statement retains the supplement\'s stated scan/certification scope.'
    return result


def frontiers():
    records=[]
    def add(identifier, group, charge, claims_ids, evidence, **extra):
        records.append({'id':identifier,'group':group,'native_input':charge,
                        'claim_ids':claims_ids,'input_references':evidence,**extra})
    for num in [5,6,7,8]:
        add('F4-power%d'%num,'F4',power(num),
            {5:['f4-power5-formal-survival'],6:['f4-power6-full-second-zero-report','f4-power6-formal-survival-report'],
             7:['f4-power7-first-zero','f4-power7-later-open'],8:['f4-power8-first-hit']}[num],
            [reference('algorithms/dual_singlets_20261118/f4/j36_source/focused_powers/sixth_survival.md','c778cd7c3ceb72ad8082d791b2aee5f0df082164')])
    for group,num,ids in [('SO(11)',6,['so11-power6-d3','so11-power6-d3-codex']),
                         ('SO(12)',7,['so12-power7-d3']),('SO(13)',7,['so13-power7-d3']),
                         ('SO(13)',8,['so13-power8-d3-product']),('SO(14)',8,['so14-power8-d3-report']),
                         ('SO(19)',10,['so19-finite-page-report']),('SO(25)',13,['so25-finite-page-report']),
                         ('SO(26)',14,['so26-finite-page-report'])]:
        add(group+'-power%d'%num,group,power(num),ids,[reference('RESEARCH_DISCUSSION.md')])
    for group,num in [('SO(11)',5),('SO(13)',6)]:
        add(group+'-one-scalar',group,native([1,0,0],[num+1,num],2*num+2),
            [group.lower().replace('(','').replace(')','')+'-one-scalar-earlier-pages'],
            [reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/one_scalar/README.md','7deb71acbdf828bd46edad5ddf503c60cbdfd58b')])
    add('SO13-quadratic-level46','SO(13)',native([2,0,0],[7,7],16),['so13-level46-two-classes'],
        [reference('algorithms/dual_singlets_20261118/orthogonal_glueballs/quadratic_so13/README.md','6a6f5a708ab1f73f00f8baafd9a133e824fe8238')],
        required_independent_classes=2)
    for length in [11,13]:
        add('G2-level36-length%d'%length,'G2',native([6,6,6],[0,0],length),
            ['g2-level36-length%d-open'%length],[reference('RESEARCH_DISCUSSION.md')])
    add('G2-power3','G2',power(3),['g2-power3-formal-report','g2-current-products-conditional'],[reference('RESEARCH_DISCUSSION.md')])
    for identifier,charge,ids in [('SU3-R007',native([7,4,4],[0,0],10,charge_conjugation='odd'),['su3-r007-formal-survival']),
        ('SU3-F1-second-seed',native([8,4,4],[0,0],11,charge_conjugation='even'),['su3-f1-formal-tower']),
        ('SU3-nine-letter',native([5,2,2],[2,2],9),['su3-nine-letter-formal-survival']),
        ('SU3-seven-letter',native([3,0,0],[4,4],7,charge_conjugation='odd',dotted_spin=0),
         ['su3-seven-letter-generator-count','su3-seven-letter-generator-explicit'])]:
        add(identifier,'SU(3)',charge,ids,[reference('RESEARCH_DISCUSSION.md')])
    return records
