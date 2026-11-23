const fs=require('fs'),path=require('path');
const {bare}=require(path.resolve(process.argv[2],'score.js'));
const {hp}=require(path.resolve(process.argv[2],'matching_score.js'));
const add=(p,w,c)=>{let d=(p.get(w)||0n)+c;if(d)p.set(w,d);else p.delete(w);};
let raw=new Map;
for(const [name,sg] of [[process.argv[3],1n],[process.argv[4],-1n]]){
 for(const line of fs.readFileSync(name,'utf8').trim().split('\n')){
  let [w,v]=line.split(/\s+/);if(w==='!')continue;
  let [s,k]=bare(w.slice(0,-1).split('.'));if(s)add(raw,k,sg*BigInt(s)*BigInt(v));
 }
}
const [projected,denominator]=hp(raw);
process.stdout.write(JSON.stringify({denominator:String(denominator),raw_terms:raw.size,projected:[...projected].sort().map(([w,c])=>[w,String(c)])})+'\n');
