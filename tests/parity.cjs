'use strict';
const fs = require('node:fs'), assert = require('node:assert/strict');
const A = require('../src/oncoflow_atlas/web/analysis.js');
let numbers = 0, maxAbsoluteError = 0;
function compare(a,b,path='$') {
  if(typeof a==='number' && typeof b==='number') {
    assert.ok(Number.isFinite(a) && Number.isFinite(b),path+' nonfinite');
    const error = Math.abs(a-b);numbers++;maxAbsoluteError=Math.max(maxAbsoluteError,error);
    assert.ok(error <= 1e-9*Math.max(1,Math.abs(b)),`${path}: ${a} != ${b}`);return;
  }
  if(a===null || b===null || typeof a!=='object' || typeof b!=='object') {assert.deepEqual(a,b,path);return;}
  assert.equal(Array.isArray(a),Array.isArray(b),path);
  assert.deepEqual(Object.keys(a).sort(),Object.keys(b).sort(),path+' keys');
  for(const key of Object.keys(b))compare(a[key],b[key],path+'.'+key);
}
const fixture=JSON.parse(fs.readFileSync(process.argv[2],'utf8')).flatMap(item=>item.views?item.views.map(v=>({...v,report:item.report})):[item]);
for(const {report,options,expected} of fixture){A.validateReport(report);compare(A.computeView(report,options),expected);}
assert.throws(()=>A.parameters({min_vaf:NaN}));
assert.throws(()=>A.parameters({top_n:3.5}));
assert.equal(A.csv([{label:'=SUM(A1:A2)',value:-2}]),"label,value\r\n'=SUM(A1:A2),-2\r\n");
assert.equal(A.escape('<img src=x onerror=alert(1)>'),'&lt;img src=x onerror=alert(1)&gt;');
console.log(JSON.stringify({status:'passed',scenarios:fixture.length,numeric_values:numbers,max_absolute_error:maxAbsoluteError,tolerance:'1e-9 relative, 1e-9 absolute'}));
