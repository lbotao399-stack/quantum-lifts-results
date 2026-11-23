// Exact ordered-trace comparison of the published SO(11) source lifts.
// The saved Gram banks contain all words with at least three trace factors.
// Every query here has six or seven, preserved by Pi and all six returns.
// This tests chosen stable derivatives, NOT the quotient by Q0(I).
const fs = require('fs');
const path = require('path');
const folder = path.resolve(process.argv[2]);
const rawFolder = process.argv[3] ? path.resolve(process.argv[3]) : folder;
const {bare,p} = require(path.join(folder,'score.js'));
const {hp} = require(path.join(folder,'matching_score.js'));
const queries = [
  ['ABC','A+','A-','B+','B-','C+','C-'],
  ['ABC','ABC','ABC','+-','+-','+-'],
];
const add=(m,k,c)=>{const v=(m.get(k)||0n)+c;if(v)m.set(k,v);else m.delete(k);};
const partition=w=>w.split('.').map(t=>t.length).sort((a,b)=>a-b).join(',');
function gcd(a,b){a=a<0n?-a:a;b=b<0n?-b:b;while(b){[a,b]=[b,a%b];}return a;}
function fraction(a,b){const d=gcd(a,b);a/=d;b/=d;return b===1n?String(a):String(a)+'/'+String(b);}
const wanted=new Set(queries.map(q=>q.map(t=>t.length).sort((a,b)=>a-b).join(',')));
const blocks=new Map([...wanted].map(k=>[k,new Map]));
const input=[];
for(const [name,coefficient] of [['cycle.raw',3n],['path3.raw',-15n],['path3-alt.raw',500n]]){
  let count=0;
  for(const line of fs.readFileSync(path.join(rawFolder,name),'utf8').trim().split('\n')){
    const [word,value]=line.split(/\s+/);if(word==='!')continue;
    const [sign,key]=bare(word.slice(0,-1).split('.'));
    if(!sign)continue;
    const block=blocks.get(partition(key));
    if(block){add(block,key,coefficient*BigInt(sign)*BigInt(value));count++;}
  }
  input.push({name,selected_raw_terms:count});
}
const output=[];
for(const [shape,raw] of blocks){
  const [projected,denominator]=hp(raw);
  let current=projected;
  const sizes=[current.size];
  for(const scalar of 'ABC')for(const ghost of ['-','+']){
    current=p(current,scalar,ghost,false);sizes.push(current.size);
  }
  for(const query of queries){
    if(query.map(t=>t.length).sort((a,b)=>a-b).join(',')!==shape)continue;
    const [sign,key]=bare(query);
    const value=sign?(current.get(key)||0n)*BigInt(sign):0n;
    output.push({query,partition:shape,raw_canonical_terms:raw.size,
      projected_terms:projected.size,return_supports:sizes,
      projector_denominator:String(denominator),
      signed_dot_derivative:fraction(value,3n*denominator)});
  }
}
process.stdout.write(JSON.stringify({scope:'Two exact coefficients of the chosen stable derivative; not a Q0I quotient test',input,output},null,2)+'\n');
