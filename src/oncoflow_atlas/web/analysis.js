/* Portable view calculations: Python analysis.py is checked against this implementation. */
(function (root) {
  'use strict';
  const defaults = {min_vaf: .05, top_n: 12, min_purity: .4, subtype: 'All', search: ''};
  const lexical = (a, b) => a < b ? -1 : a > b ? 1 : 0;
  const sum = values => values.reduce((a, b) => a + b, 0);
  const mean = values => values.length ? sum(values) / values.length : 0;
  const variance = values => {if(values.every(v=>v===values[0]))return 0;const m = mean(values); return values.length > 1 ? sum(values.map(v => (v - m) ** 2)) / (values.length - 1) : 0;};
  const median = values => {if (!values.length) return null; const a = [...values].sort((x,y) => x-y), i = Math.floor(a.length/2); return a.length % 2 ? a[i] : (a[i-1]+a[i])/2;};
  const count = (values, key) => {const m = new Map(); values.forEach(v => m.set(v[key], (m.get(v[key]) || 0)+1)); return m;};
  const escape = value => String(value ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  function parameters(options = {}) {
    const p = {...defaults, ...options};
    if (Object.keys(p).some(k => !(k in defaults))) throw Error('Unknown analysis parameter');
    for (const k of ['min_vaf','min_purity']) if (typeof p[k] !== 'number' || !Number.isFinite(p[k]) || p[k] < 0 || p[k] > 1) throw Error('Invalid '+k);
    if (!Number.isInteger(p.top_n) || p.top_n < 3 || p.top_n > 50) throw Error('Invalid top_n');
    if (typeof p.subtype !== 'string' || typeof p.search !== 'string') throw Error('Invalid cohort filter');
    p.search = p.search.trim().toLowerCase(); return p;
  }
  function validateReport(r) {
    const check = (test, message) => {if (!test) throw Error(message);};
    const text = value => typeof value === 'string' && value.trim().length > 0;
    const finite = value => typeof value === 'number' && Number.isFinite(value);
    check(r && r.schema_version === '2.0' && text(r.cohort), 'Expected a schema-2.0 analysis report');
    check(Array.isArray(r.metadata) && r.metadata.length >= 3 && r.metadata.length <= 500, 'Invalid sample count');
    const ids = r.metadata.map(m => m.sample_id), idset = new Set(ids);
    check(ids.every(text) && idset.size === ids.length, 'Invalid or duplicate sample identifiers');
    check(new Set(r.metadata.map(m => m.patient_id)).size === ids.length, 'Patient identifiers must be distinct');
    r.metadata.forEach(m => check(text(m.subtype) && text(m.patient_id) && (m.purity === null || (finite(m.purity) && m.purity >= 0 && m.purity <= 1)), 'Invalid metadata or purity'));
    function matrix(m, name, discrete) {
      check(m && Array.isArray(m.genes) && m.genes.length >= (discrete ? 1 : 2) && m.genes.length <= 60000 && m.genes.every(text) && new Set(m.genes).size === m.genes.length, 'Invalid '+name+' genes');
      check(Array.isArray(m.samples) && m.samples.length === ids.length && m.samples.every((s,i) => s === ids[i]), 'Misaligned '+name+' samples');
      check(Array.isArray(m.values) && m.values.length === m.genes.length && m.genes.length*ids.length <= 20000000, 'Invalid '+name+' dimensions');
      m.values.forEach(row => check(Array.isArray(row) && row.length === ids.length && row.every(v => finite(v) && (discrete ? [-2,-1,0,1,2].includes(v) : v >= 0)), 'Invalid '+name+' values'));
    }
    matrix(r.expression, 'expression', false); check(ids.every((s,i) => r.expression.values.some(row => row[i] > 0)), 'Zero expression sample'); if (r.cna !== null) matrix(r.cna, 'CNA', true);
    check(Array.isArray(r.variants) && r.variants.length <= 1000000, 'Invalid variant records');
    check(new Set(r.variants.map(v => v.variant_id)).size === r.variants.length, 'Duplicate variant identifiers');
    r.variants.forEach(v => check(text(v.variant_id) && text(v.gene) && text(v.impact) && text(v.variant_class) && idset.has(v.sample_id) && (v.vaf === null || (finite(v.vaf) && v.vaf >= 0 && v.vaf <= 1)), 'Invalid variant or VAF'));
    check(r.analysis && Array.isArray(r.analysis.pca) && r.analysis.pca.length === ids.length, 'Invalid PCA coverage');
    check(new Set(r.analysis.pca.map(p => p.sample_id)).size === ids.length, 'Duplicate PCA sample');
    r.analysis.pca.forEach(p => check(idset.has(p.sample_id) && finite(p.PC1) && finite(p.PC2), 'Invalid PCA coordinates'));
    const ev = r.analysis.pca_variance_explained_pct;
    check(ev && finite(ev.PC1) && finite(ev.PC2) && ev.PC1 >= 0 && ev.PC2 >= 0 && ev.PC1+ev.PC2 <= 100.000001, 'Invalid PCA variance');
    check(Array.isArray(r.analysis.correlation) && r.analysis.correlation.length === ids.length && r.analysis.correlation.every(row => Array.isArray(row) && row.length === ids.length && row.every(v => v === null || (finite(v) && Math.abs(v) <= 1.000001))), 'Invalid correlation matrix');
    check(r.source && typeof r.source === 'object' && r.inputs && typeof r.inputs === 'object' && /^[a-f0-9]{64}$/.test(r.data_fingerprint), 'Missing provenance');
    parameters(r.parameters); return r;
  }
  function computeView(r, overrides = {}) {
    const p = parameters({...r.parameters, ...overrides});
    const meta = r.metadata.filter(m => (p.subtype === 'All' || m.subtype === p.subtype) && (!p.search || m.sample_id.toLowerCase().includes(p.search)));
    const ids = meta.map(m => m.sample_id), set = new Set(ids), indices = ids.map(s => r.expression.samples.indexOf(s));
    const records = r.variants.filter(v => set.has(v.sample_id));
    const kept = records.filter(v => v.vaf !== null && v.vaf >= p.min_vaf);
    const rows = r.expression.values.map(row => indices.map(i => Math.log2(row[i]+1)));
    const variances = rows.map(variance);
    const order = ids.length ? r.expression.genes.map((g,i) => i).sort((a,b) => variances[b]-variances[a] || lexical(r.expression.genes[a], r.expression.genes[b])).slice(0,p.top_n) : [];
    const top = order.map(i => ({gene:r.expression.genes[i], variance:variances[i]}));
    const heat = order.map(i => {const m = mean(rows[i]), sd = Math.sqrt(variances[i]); return {gene:r.expression.genes[i], z:rows[i].map(v => sd > 0 ? (v-m)/sd : 0)};});
    const geneCounts = count(kept,'gene'), sampleCounts = count(kept,'sample_id'), missingVaf = count(records.filter(v => v.vaf === null),'sample_id'), altered = new Map();
    kept.forEach(v => {if (!altered.has(v.gene)) altered.set(v.gene,new Set()); altered.get(v.gene).add(v.sample_id);});
    const varSummary = [...geneCounts].sort((a,b) => b[1]-a[1] || lexical(a[0],b[0])).map(([g,n]) => ({gene:g,records:n,samples:altered.get(g).size}));
    const cnaRows = r.cna ? r.cna.values.map(row => indices.map(i => row[i])) : null;
    const cnaSummary = r.cna && ids.length ? r.cna.genes.map((g,i) => ({gene:g,amplified_samples:cnaRows[i].filter(v => v === 2).length,deep_deleted_samples:cnaRows[i].filter(v => v === -2).length})).filter(row => row.amplified_samples+row.deep_deleted_samples > 0).sort((a,b) => (b.amplified_samples+b.deep_deleted_samples)-(a.amplified_samples+a.deep_deleted_samples) || lexical(a.gene,b.gene)) : [];
    const metrics = meta.map((m,i) => {
      const reasons = [];
      if (m.purity === null) reasons.push('Purity not reported'); else if (m.purity < p.min_purity) reasons.push('Purity below threshold');
      if (missingVaf.get(m.sample_id)) reasons.push('Variant VAF unavailable');
      return {...m,variant_count:sampleCounts.get(m.sample_id)||0,variants_missing_vaf:missingVaf.get(m.sample_id)||0,
        expression_signal:mean(rows.map(row => row[i])),detected_genes:r.expression.values.filter(row => row[indices[i]] > 0).length,
        amplified_genes:cnaRows ? cnaRows.filter(row => row[i] === 2).length : null,
        deep_deleted_genes:cnaRows ? cnaRows.filter(row => row[i] === -2).length : null,
        qc_status:reasons.length ? 'REVIEW' : 'PASS',qc_reasons:reasons};
    });
    const purities = meta.map(m => m.purity).filter(v => v !== null);
    const subtypes = [...new Set(meta.map(m => m.subtype))].sort(lexical);
    const groups = subtypes.map(type => {const group = metrics.filter(m => m.subtype === type);return {subtype:type,samples:group.length,median_purity:median(group.map(m => m.purity).filter(v => v !== null)),median_variants:median(group.map(m => m.variant_count)),mean_expression_signal:mean(group.map(m => m.expression_signal))};});
    const unknown = records.filter(v => v.vaf === null).length, lowPurity = purities.filter(v => v < p.min_purity).length;
    const qc = [
      {metric:'Input identifiers and sample alignment',value:'Validated',status:'PASS'},
      {metric:'Missing expression cells',value:0,status:'PASS'},
      {metric:'Samples without reported purity',value:meta.length-purities.length,status:meta.length > purities.length ? 'NOT_REPORTED' : 'PASS'},
      {metric:'Known purity below threshold',value:lowPurity,status:lowPurity ? 'REVIEW' : 'PASS'},
      {metric:'Records with unavailable VAF',value:unknown,status:unknown ? 'REVIEW' : 'PASS'},
      {metric:'Known VAF below threshold',value:records.filter(v => v.vaf !== null && v.vaf < p.min_vaf).length,status:'FILTERED'},
      {metric:'Copy-number profile',value:r.cna ? r.cna.genes.length : 'Not supplied',status:r.cna ? 'PASS' : 'NOT_REPORTED'}
    ];
    const assessable = qc.filter(q => ['PASS','REVIEW'].includes(q.status));
    const summary = {samples:ids.length,patients:ids.length,genes:r.expression.genes.length,variants_input:records.length,variants_retained:kept.length,
      variant_genes:geneCounts.size,median_purity:median(purities),purity_reported:purities.length,samples_flagged:metrics.filter(m => m.qc_status === 'REVIEW').length,
      subtypes:Object.fromEntries(subtypes.map(t => [t,meta.filter(m => m.subtype === t).length])),qc_pass_rate_pct:100*assessable.filter(q => q.status === 'PASS').length/assessable.length,
      qc_rate_scope:'Assessable cohort checks only; unavailable fields are listed separately',parameters:p,pca_variance_explained_pct:r.analysis.pca_variance_explained_pct};
    return {parameters:p,sample_ids:ids,summary,sample_metrics:metrics,top_genes:top,heatmap:heat,variant_summary:varSummary,cna_summary:cnaSummary,subtype_summary:groups,qc,retained_variant_ids:kept.map(v => v.variant_id)};
  }
  function csv(rows) {
    if (!rows.length) return '';
    const keys = Object.keys(rows[0]);
    function cell(v) {
      if (Array.isArray(v)) v = v.join('; ');
      if (v === null || v === undefined) v = '';
      if (typeof v === 'string' && /^[=+@\-\t\r]/.test(v)) v = "'"+v;
      const text = String(v); return /[",\r\n]/.test(text) ? '"'+text.replace(/"/g,'""')+'"' : text;
    }
    return [keys.map(cell).join(','),...rows.map(row => keys.map(k => cell(row[k])).join(','))].join('\r\n')+'\r\n';
  }
  const api = {computeView, validateReport, parameters, escape, csv, median, variance};
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  root.OncoAnalysis = api;
})(typeof globalThis !== 'undefined' ? globalThis : this);
