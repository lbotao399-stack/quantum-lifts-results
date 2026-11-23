const fs=require('fs'),path=require('path'),assert=require('assert');
const {bare}=require(path.resolve(process.argv[2],'score.js'));
const add=(p,w,c)=>{let d=(p.get(w)||0n)+c;if(d)p.set(w,d);else p.delete(w);};
const pair={a:'BC',b:'AC',c:'AB'};
function replace(c,h,k){
 if(c===k)return [1,h];
 if(!pair[c]||!pair[c].includes(k))return [0,''];
 let z=pair[c].replace(k,h);if(z[0]===z[1])return [0,''];
 let s=z[0]<z[1]?1:-1;z=[...z].sort().join('');
 return [s,Object.keys(pair).find(x=>pair[x]===z)];
}
function ef(P,h,k){
 let Q=new Map;
 for(let [w,c]of P){let words=w.split('.');
  for(let i=0;i<words.length;i++)for(let j=0;j<words[i].length;j++){
   let [s,v]=replace(words[i][j],h,k);if(!s)continue;
   let changed=[...words];changed[i]=words[i].slice(0,j)+v+words[i].slice(j+1);
   let [a,b]=bare(changed);if(a)add(Q,b,BigInt(a*s)*c);
  }
 }return Q;
}
function ext(P,h,k,r){
 let den=1n,F=[1n];for(let j=1;j<=5;j++){den*=BigInt(j*(r+j));F.push(den);}
 let Q=new Map([...P].map(([w,c])=>[w,c*den])),R=P;
 for(let j=1;j<=5;j++){R=ef(R,h,k);if(!R.size)break;let V=R;for(let t=0;t<j;t++)V=ef(V,k,h);
  let z=(j%2?-den:den)/F[j];for(let [w,c]of V)add(Q,w,z*c);
 }return [Q,den];
}
function project(P){let D=1n;for(let [h,k,r]of [['B','C',1],['A','C',2],['A','B',1]]){let [Q,d]=ext(P,h,k,r);P=Q;D*=d;}
 assert.equal(ef(P,'B','C').size,0);assert.equal(ef(P,'A','B').size,0);return [P,D];}
let raw=new Map;
for(let line of fs.readFileSync(process.argv[3],'utf8').trim().split('\n')){
 let [w,c]=line.split(/\s+/);if(w==='!')continue;
 let [s,v]=bare(w.replaceAll('D','a').slice(0,-1).split('.'));if(s)add(raw,v,BigInt(s)*BigInt(c));
}
let [projected,D]=project(raw);
let image=new Map;
for(let [w,c] of projected){let words=w.split('.');
 for(let i=0;i<words.length;i++)for(let j=0;j<words[i].length;j++){
  let v=pair[words[i][j]];if(!v)continue;
  for(let [p,s]of [[v,1],[[...v].reverse().join(''),-1]]){
   let changed=[...words];changed[i]=words[i].slice(0,j)+p+words[i].slice(j+1);
   let [a,b]=bare(changed);if(a)add(image,b,c*BigInt(a*s));
  }
 }
}
process.stdout.write(JSON.stringify({denominator:String(D),raw_terms:raw.size,
 projected:[...projected].sort().map(([w,c])=>[w,String(c)]),
 Q0_projected:[...image].sort().map(([w,c])=>[w,String(c)])})+'\n');
