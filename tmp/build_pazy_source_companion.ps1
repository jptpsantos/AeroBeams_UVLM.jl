param(
    [string]$SnapshotDate = [datetime]::Now.ToString('d MMMM yyyy', [cultureinfo]::InvariantCulture)
)

$ErrorActionPreference = 'Stop'
$repo = (Get-Location).Path
$utf8 = New-Object System.Text.UTF8Encoding($false)
$exampleRoot = 'lib/WingPropellerUVLM/examples/pazy_aerobeams_coupled/'
$entryNames = @(
    'run_pazy_wing_uvlm.jl', 'PazyWingUVLMCoupling.jl',
    'PazyWingUVLMVisualization.jl', 'run_pazy_stability_sweep.jl',
    'PazyStabilitySweepTools.jl', 'plot_pazy_flutter_boundaries.jl',
    'run_pazy_wing_uvlm_static.jl', 'PazyWingUVLMStatic.jl',
    'run_aero_comparison.jl'
)
$paths = [System.Collections.Generic.List[string]]::new()
foreach ($name in $entryNames) { $paths.Add($exampleRoot + $name) }
$paths.Add('lib/WingPropellerUVLM/studies/src/MovingBlockDamping.jl')
foreach ($directory in @('src', 'lib/WingPropellerUVLM/src')) {
    foreach ($item in (Get-ChildItem -LiteralPath (Join-Path $repo $directory) -Recurse -File -Filter '*.jl' | Sort-Object FullName)) {
        $paths.Add($item.FullName.Substring($repo.Length + 1).Replace('\','/'))
    }
}
$records = [System.Collections.Generic.List[object]]::new()
$declarationCount = 0
$lineCount = 0
foreach ($relative in $paths) {
    $full = Join-Path $repo $relative
    $source = [System.IO.File]::ReadAllText($full, [System.Text.Encoding]::UTF8)
    $lines = $source -split '\r?\n'
    $functions = [System.Collections.Generic.List[object]]::new()
    for ($i = 0; $i -lt $lines.Length; $i++) {
        $name = $null
        if ($lines[$i] -match '^\s*function\s+([^\s(]+)') {
            $name = $Matches[1]
        } elseif ($lines[$i] -match '^\s*([\p{L}_][\p{L}\p{N}_.!∂]*)\([^=]*\)\s*(?:where\s+[^=]+)?\s*=') {
            $name = $Matches[1]
        }
        if ($null -ne $name) {
            $functions.Add([ordered]@{name=$name; line=$i+1})
            $declarationCount++
        }
    }
    $group = if ($relative.StartsWith($exampleRoot)) { 'Pazy examples' }
             elseif ($relative.Contains('MovingBlockDamping')) { 'Damping analysis' }
             elseif ($relative.StartsWith('src/')) { 'AeroBeams package' }
             else { 'WingPropellerUVLM package' }
    $records.Add([ordered]@{
        path=$relative; group=$group; source=$source;
        functions=@($functions.ToArray()); lines=$lines.Length;
        sha256=(Get-FileHash -LiteralPath $full -Algorithm SHA256).Hash.ToLowerInvariant()
    })
    $lineCount += $lines.Length
}
$json = ConvertTo-Json -InputObject @($records.ToArray()) -Depth 6 -Compress
# Keep embedded source inert even if a Julia string contains HTML markup.
$json = $json.Replace('<','\u003c').Replace('>','\u003e').Replace('&','\u0026')
$template = @'
<!doctype html>
<html lang="en">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pazy UVLM + AeroBeams — source companion</title>
<style>
:root{color-scheme:light;--paper:#fbfcfd;--ink:#142235;--muted:#53647a;--line:#d5dfe8;--blue:#175ca5;--accent:#e8f1fb}
*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.55 system-ui,Segoe UI,sans-serif}
header{padding:22px 28px;border-bottom:1px solid var(--line);background:#fff}h1{margin:0 0 6px;font-size:24px;letter-spacing:-.5px}p{margin:6px 0;color:var(--muted)}a{color:var(--blue)}
.layout{display:grid;grid-template-columns:360px minmax(0,1fr);height:calc(100vh - 184px);min-height:500px}
aside{padding:18px;border-right:1px solid var(--line);overflow:auto;background:#f1f5f9}input{width:100%;padding:11px 12px;border:1px solid #adbdcd;border-radius:7px;font:inherit;background:white}
#counts{font-size:12px;margin:8px 0 15px}h2{font-size:13px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:20px 0 8px}
button{font:inherit;color:var(--blue);background:transparent;border:0;text-align:left;cursor:pointer;padding:4px 0;overflow-wrap:anywhere}button:hover{text-decoration:underline}
.file{margin:7px 0 12px}.file>button{font-size:14px;font-weight:600}.path{font-size:11px;color:var(--muted);overflow-wrap:anywhere}.func{display:block;font:12px/1.5 ui-monospace,Consolas,monospace;padding:2px 7px}.func:before{content:'↳ ';color:#6e8298}
details{margin-top:4px}summary{font-size:12px;color:var(--muted);cursor:pointer}main{display:flex;flex-direction:column;min-width:0;background:#fff}.toolbar{padding:14px 20px;border-bottom:1px solid var(--line)}#filename{font:600 15px/1.4 ui-monospace,Consolas,monospace;overflow-wrap:anywhere}#metadata{font-size:11px;overflow-wrap:anywhere}
#source{overflow:auto;flex:1;margin:0;padding:12px 0 40px;font:13px/1.6 ui-monospace,Consolas,monospace;tab-size:4;counter-reset:ln}.line{display:flex;min-width:max-content;min-height:21px;padding-right:24px}.number{display:inline-block;flex:0 0 66px;width:66px;text-align:right;padding-right:16px;user-select:none;color:#8b99a8;text-decoration:none}.code{white-space:pre}.line:target,.line.selected{background:#fff0bd}.line:hover{background:#f2f7fc}.line.selected:hover{background:#fff0bd}
@media(max-width:850px){.layout{grid-template-columns:270px minmax(0,1fr)}header{padding:16px}.layout{height:calc(100vh - 210px)}}
@media(max-width:600px){.layout{display:block;height:auto}aside{max-height:300px}main{height:75vh}h1{font-size:20px}}
</style>
<header>
  <h1>Pazy UVLM + AeroBeams · Source companion</h1>
  <p>Complete local Julia source captured on __SNAPSHOT_DATE__. <a href="pazy-uvlm-aerobeams-walkthrough.md">Read the detailed walkthrough</a>.</p>
  <p>This is a source snapshot, including working-copy edits. Package source includes available and conditional routines; the walkthrough identifies the Pazy execution path. The declaration index is a text-based navigation aid, not a measured call graph. Third-party package implementations are not included.</p>
</header>
<div class="layout">
 <aside><input id="search" type="search" placeholder="Find a file, function or source text…" aria-label="Find source"><p id="counts"></p><nav id="navigation" aria-label="Source files and functions"></nav></aside>
 <main><div class="toolbar"><div id="filename"></div><p id="metadata"></p></div><pre id="source" aria-label="Numbered source code"></pre></main>
</div>
<script type="application/json" id="data">__SOURCE_DATA__</script>
<script>
'use strict';
const files=JSON.parse(document.getElementById('data').textContent);
const nav=document.getElementById('navigation'), query=document.getElementById('search'), source=document.getElementById('source');
let selected=-1;
function openSource(index,line=1,push=true){
 const f=files[index]; if(!f)return;
 if(selected!==index){
   document.getElementById('filename').textContent=f.path;
   document.getElementById('metadata').textContent=`${f.group} · ${f.lines.toLocaleString()} lines · SHA-256 ${f.sha256}`;
   const fragment=document.createDocumentFragment();
   f.source.split(/\r?\n/).forEach((text,i)=>{
      const row=document.createElement('span');row.className='line';row.id='L'+(i+1);
      const num=document.createElement('a');num.className='number';num.textContent=i+1;num.href='#';
      num.addEventListener('click',e=>{e.preventDefault();openSource(index,i+1)});
      const code=document.createElement('span');code.className='code';code.textContent=text||' ';
      row.append(num,code);fragment.append(row);
   });source.replaceChildren(fragment);selected=index;
 }
 document.querySelectorAll('.line.selected').forEach(r=>r.classList.remove('selected'));
 const row=document.getElementById('L'+line);if(row){row.classList.add('selected');row.scrollIntoView({block:'center'});}
 if(push)history.replaceState(null,'','#file='+encodeURIComponent(f.path)+'&line='+line);
}
function renderNavigation(){
 const q=query.value.toLocaleLowerCase().trim();nav.replaceChildren();let count=0,currentGroup='';
 files.forEach((f,index)=>{
   const pathMatch=f.path.toLocaleLowerCase().includes(q);
   const funcs=f.functions.filter(fn=>fn.name.toLocaleLowerCase().includes(q));
   const bodyMatch=q && f.source.toLocaleLowerCase().includes(q);
   if(q&&!pathMatch&&!funcs.length&&!bodyMatch)return;
   if(currentGroup!==f.group){currentGroup=f.group;const h=document.createElement('h2');h.textContent=f.group;nav.append(h);}
   count++;const section=document.createElement('div');section.className='file';
   const button=document.createElement('button');button.textContent=f.path.split('/').pop();
   const firstLine=q&&!pathMatch&&!funcs.length?f.source.split(/\r?\n/).findIndex(l=>l.toLocaleLowerCase().includes(q))+1:1;
   button.addEventListener('click',()=>openSource(index,Math.max(1,firstLine)));section.append(button);
   const p=document.createElement('div');p.className='path';p.textContent=f.path;section.append(p);
   const list=q&&!pathMatch?funcs:f.functions;
   if(list.length){const d=document.createElement('details');d.open=!!q;const s=document.createElement('summary');s.textContent=list.length+' declarations';d.append(s);
      list.forEach(fn=>{const b=document.createElement('button');b.className='func';b.textContent=fn.name+' · '+fn.line;b.addEventListener('click',()=>openSource(index,fn.line));d.append(b)});section.append(d);}
   nav.append(section);
 });
 document.getElementById('counts').textContent=`${count} / ${files.length} files · __TOTAL_LINES__ source lines · __TOTAL_DECLARATIONS__ indexed declarations`;
}
function fromHash(){const p=new URLSearchParams(location.hash.slice(1));const index=files.findIndex(f=>f.path===p.get('file'));openSource(index<0?0:index,Math.max(1,Number(p.get('line'))||1),false)}
query.addEventListener('input',renderNavigation);window.addEventListener('hashchange',fromHash);renderNavigation();fromHash();
</script>
</html>
'@
$html = $template.Replace('__SOURCE_DATA__',$json).Replace('__SNAPSHOT_DATE__', [System.Net.WebUtility]::HtmlEncode($SnapshotDate)).Replace('__TOTAL_LINES__', $lineCount.ToString('N0',[cultureinfo]::InvariantCulture)).Replace('__TOTAL_DECLARATIONS__',$declarationCount.ToString('N0',[cultureinfo]::InvariantCulture))
$output = Join-Path $repo 'docs/src/pazy-uvlm-aerobeams-source.html'
[System.IO.File]::WriteAllText($output, $html, $utf8)
[pscustomobject]@{Path=$output; Files=$records.Count; Lines=$lineCount; Declarations=$declarationCount; Bytes=(Get-Item -LiteralPath $output).Length}
