'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const {pathToFileURL}=require('node:url');
let playwright;try{playwright=require('playwright');}catch{playwright=require(process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES+'/playwright');}
const root=path.resolve(__dirname,'..'),proof=path.join(root,'validation'),assets=path.join(root,'assets');
const capture=process.argv.includes('--capture'),checks=[],errors=[],requests=[];
const sha=filename=>crypto.createHash('sha256').update(fs.readFileSync(filename)).digest('hex');
function passed(name){checks.push(name);}
async function screenshot(page,name,fullPage=false){const file=path.join(assets,name);fs.writeFileSync(file+'.tmp',await page.screenshot({fullPage}));fs.renameSync(file+'.tmp',file);}
async function navigate(page,route){await page.locator(`[data-route="${route}"]`).click();await page.waitForFunction(route=>document.querySelector(`[data-route="${route}"]`).getAttribute('aria-current')==='page',route);}
(async()=>{
  fs.mkdirSync(proof,{recursive:true});fs.mkdirSync(assets,{recursive:true});
  const launch={headless:true,args:['--no-sandbox']};
  if(process.env.ONCOFLOW_BROWSER)launch.executablePath=process.env.ONCOFLOW_BROWSER;
  const browser=await playwright.chromium.launch(launch);
  try{
    const context=await browser.newContext({viewport:{width:1440,height:1000},offline:true});
    const page=await context.newPage();
    page.on('pageerror',e=>errors.push(e.message));
    page.on('console',m=>{if(m.type()==='error')errors.push(m.text());});
    page.on('request',r=>{if(/^https?:/.test(r.url()))requests.push(r.url());});
    await page.goto(pathToFileURL(path.join(root,'docs/index.html')).href);
    await page.waitForFunction(()=>document.querySelector('#view-status').textContent.startsWith('36 samples'));
    assert.match(await page.locator('#view-status').innerText(),/3,592 retained/);passed('Public dataset starts with 36 samples and 3,592 retained records');
    assert.equal(await page.locator('#dataset option').count(),2);passed('Public and synthetic reports are separate datasets');
    const routes=['overview','expression','variants','copy-number','gene','quality','provenance'];
    for(let i=0;i<routes.length;i++){
      await navigate(page,routes[i]);
      const text=await page.locator('#workspace').innerText();assert.ok(text.length>250);assert.ok(!/\bNaN\b/.test(text),routes[i]+' contains a nonfinite display value');
      passed('Desktop route '+routes[i]);
      if(capture){await screenshot(page,`tour-${i+1}.png`);if(['overview','expression','gene'].includes(routes[i]))await screenshot(page,`workbench-${routes[i]}.png`,true);}
    }
    await navigate(page,'gene');await page.locator('#gene-input').fill('IDH1');await page.locator('#inspect-gene').click();
    assert.match(await page.locator('#workspace').innerText(),/IDH1.*RNA profile present.*CNA profile present.*5 retained/s);passed('IDH1 joins RNA, CNA and five retained records');
    await navigate(page,'overview');
    await page.locator('#vaf').evaluate(el=>{el.value='0.3';el.dispatchEvent(new Event('input',{bubbles:true}));});await page.locator('#apply').click();
    assert.match(await page.locator('#view-status').innerText(),/2,079 retained/);passed('Exact 0.30 VAF boundary retains 2,079 records');
    await page.locator('#subtype').selectOption('Perihilar Cholangiocarcinoma');await page.locator('#apply').click();
    assert.match(await page.locator('#view-status').innerText(),/^4 samples/);assert.equal(await page.locator('#pca-chart circle').count(),4);passed('Histology filter changes summary and PCA point selection');
    const downloadPromise=page.waitForEvent('download');await page.locator('#export-csv').click();const download=await downloadPromise;
    const csv=fs.readFileSync(await download.path(),'utf8').trim().split(/\r?\n/);assert.equal(csv.length,5);assert.ok(csv.slice(1).every(row=>row.includes('Perihilar Cholangiocarcinoma')));passed('CSV export contains only four selected samples');
    await page.locator('#reset').click();const sample=JSON.parse(fs.readFileSync(path.join(root,'results/public/analysis.json'))).metadata[0].sample_id;
    await page.locator('#sample-search').fill(sample);await page.locator('#apply').click();assert.match(await page.locator('#view-status').innerText(),/^1 samples/);passed('One-sample selection renders');
    await navigate(page,'expression');assert.equal(await page.locator('#expression-heatmap td').count(),12);assert.ok((await page.locator('#expression-heatmap td').evaluateAll(nodes=>nodes.map(n=>n.title))).every(t=>t.includes('z = 0.0000')));passed('One-sample heatmap has defined zero z-scores');
    await page.locator('#sample-search').fill('NO-SUCH-SAMPLE');await page.locator('#apply').click();assert.match(await page.locator('#view-status').innerText(),/^0 samples/);assert.match(await page.locator('#workspace').innerText(),/No cells/);passed('Empty selection is explicit');
    await page.locator('#reset').click();await page.locator('#dataset').selectOption('1');await navigate(page,'copy-number');assert.match(await page.locator('#workspace').innerText(),/does not supply copy-number/);passed('Synthetic dataset does not invent a CNA assay');
    await navigate(page,'quality');assert.match(await page.locator('#source-badge').innerText(),/Authored/);passed('Synthetic provenance stays visible');
    const payload=JSON.parse(fs.readFileSync(path.join(root,'results/example/analysis.json')));
    const invalid=structuredClone(payload);invalid.variants[0].sample_id='unknown';
    await page.locator('#import-file').setInputFiles({name:'invalid.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(invalid))});
    await page.waitForFunction(()=>document.querySelector('#toast').textContent.startsWith('Invalid report:'));assert.equal(await page.locator('#dataset option').count(),2);passed('Invalid import is rejected without adding a dataset');
    payload.cohort='<img src=x onerror="window.attack=true">';payload.source.name=payload.cohort;payload.source.url='javascript:window.attack=true';
    await page.locator('#import-file').setInputFiles({name:'escaped.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(payload))});
    await page.waitForFunction(()=>document.querySelector('#dataset').options.length===3);await navigate(page,'provenance');
    assert.equal(await page.locator('#workspace img').count(),0);assert.equal(await page.locator('a[href^="javascript:"]').count(),0);assert.equal(await page.evaluate(()=>window.attack),undefined);passed('Imported labels and unsafe URLs cannot execute HTML or script');
    await page.locator('#theme').click();assert.equal(await page.locator('html').getAttribute('data-theme'),'light');await page.reload();await page.waitForFunction(()=>document.querySelector('#view-status').textContent.startsWith('36 samples'));assert.equal(await page.locator('html').getAttribute('data-theme'),'light');passed('Theme choice persists');
    await page.locator('#theme').click();await navigate(page,'overview');
    await page.keyboard.press('Tab');assert.ok(await page.evaluate(()=>document.activeElement!==document.body));passed('Keyboard focus is available');
    for(const width of [390,360]){
      await page.setViewportSize({width,height:844});
      for(const route of routes){await navigate(page,route);assert.ok(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth+1),`${route} overflows at ${width}px`);passed(`Mobile ${width}px route ${route}`);}
    }
    if(capture){await page.setViewportSize({width:390,height:844});await navigate(page,'overview');await screenshot(page,'workbench-mobile.png',true);}
    assert.deepEqual(errors,[]);assert.deepEqual(requests,[]);passed('No console errors or external data requests in offline mode');
    const result={status:'passed',browser:await browser.version(),playwright:require(playwright===undefined?'playwright/package.json':(process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES?process.env.CODEX_PRIMARY_RUNTIME_NODE_MODULES+'/playwright/package.json':'playwright/package.json')).version,checks:checks.length,passed:checks,console_errors:errors,external_requests:requests,docs_sha256:sha(path.join(root,'docs/index.html'))};
    fs.writeFileSync(path.join(proof,'browser-checks.json'),JSON.stringify(result,null,2)+'\n');console.log(JSON.stringify(result,null,2));
    await context.close();
  }finally{await browser.close();}
})().catch(e=>{console.error(e.stack);process.exit(1);});
