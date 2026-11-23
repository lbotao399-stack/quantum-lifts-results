exec(open(__file__.replace('character_frontier.py','character_scan.py')).read().split('start=time.monotonic()')[0])
start=time.monotonic()
records=[]
for rank,k in [(4,7),(5,9),(4,6),(4,8)]:
 w=(3,3,k,0,0);L=k+4
 row={'rank':rank,'weight':w,'length':L,'cochains':{}}
 for g,n in [('Sp',2*rank),('SO',2*rank+1),('Sp',2*(L+1))]:
  row['cochains'][g+str(n)]={ell:dimensions(w,ell,n,g) for ell in range(L-1,L+2)}
 print(json.dumps(row),flush=True);records.append(row)
 (HERE/'character_frontier.json').write_text(json.dumps({'records':records,'seconds':time.monotonic()-start},indent=2))
