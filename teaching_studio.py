"""Interactive client-side playback for algorithm teaching traces.

Only one user-selected read is animated. Benchmark search remains server-side.
"""
import html
import json
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from search_engine import METHODS, Matcher, Target, lps_table
from teaching_trace import build_teaching_trace


EXPLAIN = {
    "Naive": "Start at every possible position. Compare left-to-right; after a mismatch slide the whole pattern by one base. Worst-case O(N×M) per query.",
    "KMP": "Build the LPS (prefix/suffix) table once. Compare left-to-right; when a partial match fails, reuse LPS and keep the DNA pointer instead of restarting. O(N+M) per query.",
    "Horspool": "Compare right-to-left. On each alignment, use the rightmost text character and the bad-character shift table to skip positions. Worst-case O(N×M).",
    "Aho-Corasick": "Build a trie from all patterns and failure links. Consume the DNA once while changing trie states; output any patterns ending at that state. O(N+Z) search after building the automaton.",
    "Python str.find": "This is a fast optimized native baseline. Python does not expose its internal comparison steps, so the animation shows completed matches only.",
    "Adaptive (ours)": "A course-level heuristic chooses an established search engine from the pattern set; this is not a new string matching algorithm.",
}


def _trace_player_html(payload):
    # text and query are validated A/C/G/T; still prevent script-tag breakouts.
    data = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    data = data.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    template = r'''<!doctype html>
<html><head><meta charset="utf-8"><style>
:root{font-family:ui-sans-serif,system-ui,-apple-system,Segoe UI,Arial,sans-serif;color:#182332;background:#fff}
*{box-sizing:border-box}body{margin:0;padding:15px}button,select,input{font:inherit}button{cursor:pointer}
#board{border:1px solid #e2e8f0;border-radius:16px;background:linear-gradient(180deg,#fff,#f8fafc);padding:16px;position:relative}
.row{display:flex;align-items:center;gap:9px;flex-wrap:wrap}.between{justify-content:space-between}.muted{font-size:12px;color:#64748b}
.title{font-size:16px;font-weight:800}.chip{font-size:11px;border-radius:99px;background:#e2e8f0;padding:4px 9px;font-weight:700}.ok{color:#166534}.danger{color:#991b1b}
.controls{background:white;border:1px solid #e2e8f0;border-radius:12px;padding:12px;margin-top:12px}
.btn{border:1px solid #cbd5e1;border-radius:9px;padding:7px 11px;background:white;color:#172033;font-weight:650;font-size:13px}
.btn.primary{background:#115e59;color:white;border-color:#115e59}.btn:disabled{opacity:.35;cursor:not-allowed}
input[type=range]{accent-color:#0f766e;min-width:170px;flex:1}
.wrap{display:grid;grid-template-columns:minmax(0,1fr) 252px;gap:13px;margin-top:13px}
.panel{background:white;border:1px solid #e2e8f0;border-radius:12px;padding:12px;min-width:0}
.smallh{font-size:12px;color:#475569;font-weight:750;text-transform:uppercase;letter-spacing:.06em;margin-bottom:8px}
.dnarow{display:grid;grid-template-columns:repeat(var(--count), minmax(0,1fr));gap:2px;min-width:550px}
.b{min-width:0;height:27px;border-radius:4px;background:#e8eef5;color:#1e293b;font:700 13px/27px ui-monospace,SFMono-Regular,Menlo,monospace;text-align:center;position:relative}
.b.before{background:#e2e8f0}.b.active{background:#facc15;color:#713f12;box-shadow:inset 0 0 0 2px #ca8a04}
.b.matched{background:#bbf7d0;color:#166534}.b.mismatch{background:#fecaca;color:#991b1b;box-shadow:inset 0 0 0 2px #ef4444}
.b.final{background:#16a34a;color:white}.b.align{outline:1px solid #64748b;outline-offset:-1px}
.b.pat{background:#dbeafe;color:#1d4ed8}.b.pat.active{background:#fb923c;color:#7c2d12}.b.blank{background:transparent;color:transparent}
.gridscroll{overflow-x:auto}.labels{display:grid;grid-template-columns:55px 1fr;gap:3px;align-items:center;margin:5px 0}
.rowtag{font-size:12px;color:#475569;font-weight:700}.indexrow .b{height:20px;line-height:20px;background:transparent;color:#64748b;font-size:9px}
.event{border-radius:10px;padding:10px 12px;background:#f1f5f9;margin-top:10px;font-size:13px;line-height:1.5;min-height:64px}
.event.good{background:#dcfce7;color:#14532d}.event.bad{background:#fee2e2;color:#7f1d1d}
.metrics{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:7px;margin-top:10px}
.metric{background:#f8fafc;border:1px solid #e2e8f0;padding:8px;border-radius:9px}
.metric strong{font-size:20px;display:block;margin-top:3px}.metric span{font-size:11px;color:#64748b}
.hititem{padding:7px;border-radius:8px;background:#f0fdf4;border:1px solid #bbf7d0;margin-bottom:7px;font-size:12px;cursor:pointer}
.hititem:hover{background:#dcfce7}.scrollhits{max-height:225px;overflow:auto}
#notice{position:absolute;right:18px;top:12px;z-index:4;display:none;max-width:310px;border-radius:13px;background:#14532d;color:white;box-shadow:0 14px 30px #0f172a40;padding:12px}
#notice .sub{color:#bbf7d0;font-size:12px}.legend{font-size:11px;color:#64748b;margin:10px 0 0}
.sw{display:inline-block;border-radius:2px;width:10px;height:10px;margin-right:4px}.overview{display:flex;flex-wrap:wrap;gap:2px;margin-top:7px;max-height:120px;overflow-y:auto;align-content:flex-start}
.base-mini{width:9px;height:9px;border-radius:2px;background:#cbd5e1;cursor:pointer}.base-mini.visited{background:#94a3b8}.base-mini.hit{background:#16a34a}.base-mini.now{outline:2px solid #f59e0b}
#modal{position:absolute;inset:0;display:none;align-items:center;justify-content:center;background:#0f172ab8;z-index:6;border-radius:14px;padding:20px}
.modalbox{background:white;border-radius:14px;padding:20px;max-width:640px;width:100%;max-height:88%;overflow:auto;box-shadow:0 20px 60px #0004}
.seqcontext{font:700 14px/2.1 ui-monospace,Menlo,monospace;word-break:break-all}.hctx{background:#bbf7d0;color:#14532d;padding:2px 0;border-radius:2px}
.kmpbox{background:#eef2ff;border-radius:8px;padding:8px;margin-top:7px;font:12px/1.5 ui-monospace,Menlo,monospace;overflow-wrap:anywhere}
@media(max-width:780px){.wrap{grid-template-columns:1fr}.dnarow{min-width:485px}.metrics{grid-template-columns:repeat(3,1fr)}#notice{top:5px;right:5px}}
</style></head><body>
<script type="application/json" id="trace-data">__DATA__</script>
<div id="board">
 <div id="notice" role="status"></div>
 <div class="row between"><div><div class="title">🔬 Live Algorithm Teaching Studio</div><div class="muted" id="subtitle"></div></div><span class="chip" id="progress"></span></div>
 <div class="controls">
  <div class="row"><button class="btn" id="prev">◀ Previous</button><button class="btn primary" id="play">▶ Play</button><button class="btn" id="next">Next ▶</button><button class="btn" id="prevhit">↤ Prev match</button><button class="btn" id="nexthit">Next match ↦</button></div>
  <div class="row" style="margin-top:12px"><span class="muted">Step</span><input type="range" id="seek" min="0" value="0"/><span class="muted" id="stepscount"></span><label class="muted">Speed <select id="delay"><option value="750">Slow</option><option value="350" selected>Normal</option><option value="120">Fast</option><option value="40">Very fast</option></select></label><label class="muted"><input type="checkbox" id="pauseOnHit" checked/> Pause at exact match</label></div>
 </div>
 <div class="wrap">
  <section class="panel">
   <div class="smallh">The pattern moves underneath the DNA</div><div class="muted" id="range"></div>
   <div class="gridscroll"><div id="dnaRow" class="labels"></div><div id="patternRow" class="labels"></div><div id="indexRow" class="labels indexrow"></div></div>
   <div class="legend"><span class="sw" style="background:#facc15"></span> compared base <span class="sw" style="background:#16a34a"></span> complete match <span class="sw" style="background:#fecaca"></span> mismatch <span class="sw" style="background:#dbeafe"></span> pattern <span class="sw" style="background:#e2e8f0"></span> unexamined</div>
   <div id="explain" class="event"></div>
   <div id="kmpinfo" class="kmpbox" style="display:none"></div>
   <div class="metrics"><div class="metric"><span>Comparisons / transitions</span><strong id="operations">0</strong></div><div class="metric"><span>Exact matches found so far</span><strong id="found">0</strong></div><div class="metric"><span>Current text index</span><strong id="textidx">0</strong></div></div>
   <div class="smallh" style="margin-top:13px">Entire read overview — green = exact hit</div><div id="overview" class="overview"></div><div class="muted" style="margin-top:4px" id="scope"></div>
  </section>
  <aside class="panel"><div class="smallh">Matches discovered while playing</div><div class="muted" id="foundlabel"></div><div id="hitlist" class="scrollhits" style="margin-top:8px"></div><button id="openmatch" class="btn" style="margin-top:10px" disabled>Open selected match popup</button><div class="smallh" style="margin-top:17px">What the algorithm just did</div><div class="muted" id="action"></div><div class="smallh" style="margin-top:14px">Algorithm explanation</div><div class="muted" id="algorithmInfo"></div></aside>
 </div>
 <div id="modal"><div class="modalbox"><div class="row between"><b>🎯 Complete DNA pattern match</b><button class="btn" id="closemodal">✕ Close</button></div><div id="modalcontent"></div></div></div>
</div>
<script>
(() => {
 const input=JSON.parse(document.getElementById('trace-data').textContent);
 const seq=input.dna, tr=input.trace, events=tr.events, targets=input.targets;
 const $=id=>document.getElementById(id);
 const esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
 let current=0,playing=false,timer=null,lastNotice=-1,selectedHit=null;
 const hits=[];events.forEach((e,i)=>e.hits.forEach(h=>hits.push({...h,step:i})));
 const matchSteps=Array.from(new Set(hits.map(h=>h.step)));
 const currentPattern=targets[0]?.sequence||'';
 $('subtitle').textContent=`${tr.effective_algorithm} • ${seq.length.toLocaleString()} DNA bases • ${targets.length} target(s) • ${events.length.toLocaleString()} events`;
 $('algorithmInfo').textContent=input.explain;
 $('seek').max=Math.max(0,events.length-1);
 $('scope').textContent=`${hits.length} exact occurrence(s) across this selected read. Full dataset results are shown in the benchmark table above. ${tr.complete?'Full teaching trace.':'Trace limited by browser memory safeguard; increase the cap only for manageable inputs.'}`;
 const MINI=document.createDocumentFragment();for(let i=0;i<seq.length;i++){const s=document.createElement('span');s.className='base-mini';s.title=`DNA position ${i}`;s.dataset.n=i;s.onclick=()=>jumpToPosition(i);MINI.appendChild(s)}$('overview').appendChild(MINI);
 function jumpToPosition(pos){let best=0,bestDiff=Infinity;events.forEach((e,j)=>{if(e.text_index===null)return;const d=Math.abs(e.text_index-pos);if(d<bestDiff){bestDiff=d;best=j;}});seek(best)}
 function matchesVisible(step){return hits.filter(h=>h.step<=step)}
 function dnaCells(e,windowStart,windowEnd){const start=windowStart,end=windowEnd;
  let seqHTML='',patHTML='',ixHTML='';
  const pat=e.kind==='transition'||(tr.effective_algorithm==='Aho-Corasick'&&e.kind==='fallback')?e.active_prefix:currentPattern;
  let alignment=e.alignment;
  if(tr.effective_algorithm==='Aho-Corasick' && alignment==null && pat){alignment=(e.text_index??0)-pat.length+1}
  const foundInEvent=e.hits||[];
  const previous=matchesVisible(current).filter(h=>h.step<current);
  for(let pos=start;pos<end;pos++){
    let cl='b';const currentCompare=pos===e.text_index && (e.kind==='compare'||e.kind==='transition');
    const justHit=foundInEvent.some(h=>h.start<=pos&&pos<h.end);
    const priorHit=previous.some(h=>h.start<=pos&&pos<h.end);
    const inAlignment=alignment!=null&&pos>=alignment&&pos<alignment+pat.length;
    if(inAlignment)cl+=' align';if(priorHit)cl+=' matched';if(justHit)cl+=' final';
    if(currentCompare)cl+=e.matched===false?' mismatch':' active';
    seqHTML+=`<span class="${cl}" title="DNA[${pos}]=${esc(seq[pos])}">${esc(seq[pos])}</span>`;
    let char='·',pc='b blank';if(alignment!==null&&pos>=alignment&&pos<alignment+pat.length){char=pat[pos-alignment];pc='b pat';if(pos===e.text_index&&e.kind==='compare')pc+=' active'}
    patHTML+=`<span class="${pc}" title="pattern alignment">${esc(char)}</span>`;
    ixHTML+=`<span class="b">${pos%10===0?pos:'·'}</span>`;
  }
  for(const [id,label,code] of [['dnaRow','DNA',seqHTML],['patternRow',tr.effective_algorithm==='Aho-Corasick'?'Trie':'Pattern',patHTML],['indexRow','Index',ixHTML]])
    $(id).innerHTML=`<div class="rowtag">${label}</div><div class="dnarow" style="--count:${end-start}">${code}</div>`;
 }
 function toast(h){const n=$('notice');n.innerHTML=`<b>🎯 Exact match: ${esc(h.name)}</b><div class="sub">${esc(h.sequence)} • [${h.start}, ${h.end}) • strand ${esc(h.strand)}</div>`;n.style.display='block';lastNotice=current;}
 function render(){if(!events.length){$('explain').textContent='No traceable events.';return}
  const e=events[current],n=seq.length;
  const center=Math.max(0,e.text_index??e.alignment??0);const len=Math.min(52,n);const maxStart=Math.max(0,n-len);
  // Keep the DNA fixed for a chunk of time so the pattern visibly slides right.
  const start=Math.min(maxStart,Math.floor(center/44)*44),end=Math.min(n,start+len);
  dnaCells(e,start,end);$('range').textContent=`Displaying positions ${start}–${end-1} / ${n-1} (full read is in the overview below)`;
  $('progress').textContent=`Step ${current+1} / ${events.length}`;$('stepscount').textContent=`${current+1} / ${events.length}`;$('seek').value=current;
  $('operations').textContent=e.operations.toLocaleString();$('found').textContent=e.found_count.toLocaleString();$('textidx').textContent=e.text_index??'—';
  $('explain').className='event'+(e.kind==='exact_match'?' good':e.kind==='compare'&&e.matched===false?' bad':'');
  $('explain').innerHTML=`<b>${esc(e.kind.replace('_',' ').toUpperCase())}</b><div>${esc(e.message)}</div>`;
  let action=`Current pattern window begins at ${e.alignment??'—'}; pattern index ${e.pattern_index??'—'}.`;
  if(e.shift!==null)action+=` Window moves +${e.shift} DNA position(s).`;
  if(e.state!==null)action+=` Trie state=${e.state}, active prefix=${e.active_prefix||'∅'}.`;
  $('action').textContent=action;
  if(tr.effective_algorithm==='KMP'){$('kmpinfo').style.display='block';$('kmpinfo').textContent='Pattern: '+currentPattern+' | LPS: ['+tr.lps.join(', ')+'] | Watch how fallback shifts the pattern while the DNA pointer stays.'}
  else if(tr.effective_algorithm==='Horspool'){$('kmpinfo').style.display='block';$('kmpinfo').textContent='Bad-character table: '+JSON.stringify(tr.shift_table)}
  else $('kmpinfo').style.display='none';
  const visible=matchesVisible(current);$('foundlabel').textContent=`${visible.length} discovered / ${hits.length} total in this trace`;
  $('hitlist').innerHTML=visible.length?visible.map((h,i)=>`<div class="hititem" data-hit="${h.step}" title="Click to jump to exact match"><b>${esc(h.name)}</b> • ${esc(h.strand)}<div>[${h.start}, ${h.end}) • ${esc(h.sequence)}</div></div>`).join(''):'<span class="muted">No full pattern matches detected yet. Continue playing or jump to the next match.</span>';
  $('hitlist').querySelectorAll('[data-hit]').forEach(el=>el.onclick=()=>seek(Number(el.dataset.hit)));
  if(visible.length){selectedHit=visible[visible.length-1];$('openmatch').disabled=false}else{$('openmatch').disabled=true;selectedHit=null}
  const hPositions=new Set(hits.filter(h=>h.step<=current).flatMap(h=>Array.from({length:h.end-h.start},(_,j)=>h.start+j)));
  [...$('overview').children].forEach((el,i)=>{el.className='base-mini'+(i<=center?' visited':'')+(hPositions.has(i)?' hit':'')+(i===center?' now':'')});
  $('prev').disabled=current<=0;$('next').disabled=current>=events.length-1;
  if(e.hits.length && lastNotice!==current){e.hits.forEach(toast)}
 }
 function seek(i){current=Math.min(events.length-1,Math.max(0,i));if(current<0)current=0;render()}
 function stop(){playing=false;if(timer){clearTimeout(timer);timer=null}$('play').textContent='▶ Play'}
 function tick(){if(!playing)return;if(current>=events.length-1){stop();return}
  current++;render();const e=events[current];if(e.hits.length&&$('pauseOnHit').checked){stop();return}
  timer=setTimeout(tick,Number($('delay').value));}
 $('play').onclick=()=>{if(!events.length)return;if(playing){stop();return}if(current>=events.length-1)current=0;playing=true;$('play').textContent='⏸ Pause';timer=setTimeout(tick,Number($('delay').value))};
 $('prev').onclick=()=>{stop();seek(current-1)};$('next').onclick=()=>{stop();seek(current+1)};
 $('seek').oninput=ev=>{stop();seek(Number(ev.target.value))};
 function jumpHit(direction){stop();const options=matchSteps.filter(k=>direction>0?k>current:k<current);if(options.length)seek(direction>0?options[0]:options[options.length-1]);}
 $('prevhit').onclick=()=>jumpHit(-1);$('nexthit').onclick=()=>jumpHit(1);
 function modalFor(h){const l=Math.max(0,h.start-20),r=Math.min(seq.length,h.end+20);
 $('modalcontent').innerHTML=`<p><b>${esc(h.name)}</b> • strand ${esc(h.strand)} • read-level coordinates <b>[${h.start}, ${h.end})</b></p><div class="seqcontext">${esc(seq.slice(l,h.start))}<span class="hctx">${esc(seq.slice(h.start,h.end))}</span>${esc(seq.slice(h.end,r))}</div><p class="muted">Green span is the complete exact sequence. These are read offsets, not chromosome coordinates.</p>`;
 $('modal').style.display='flex';stop()}
 $('openmatch').onclick=()=>{if(selectedHit)modalFor(selectedHit)};$('closemodal').onclick=()=>{$('modal').style.display='none'};
 $('modal').onclick=e=>{if(e.target===$('modal'))$('modal').style.display='none'};
 $('prev').disabled=!events.length;$('next').disabled=!events.length;
 render();
})();
</script>
</body></html>'''
    return template.replace("__DATA__", data)


def _record_selector(reads, targets, algorithm):
    """Select real read from cached benchmark sample; no fabricated data."""
    if not reads:
        return None
    # Use one target here to keep read selection comprehensible. 10k scans are
    # inexpensive but avoid an unbounded selection menu on very large uploads.
    matcher = Matcher("Aho-Corasick", targets) if len(targets) > 1 else Matcher("KMP", targets)
    candidates = []
    for ix, r in enumerate(reads[:20000]):
        positions, _ = matcher.scan(r["seq"])
        if positions:
            candidates.append((ix, len(positions)))
            if len(candidates) >= 150:
                break
    if not candidates:
        return None
    chosen = st.selectbox(
        "Select an actual FASTQ/FASTA read containing this pattern",
        candidates,
        format_func=lambda c: f"{reads[c[0]]['source']} · {reads[c[0]]['id']} ({c[1]} matching occurrence(s), {len(reads[c[0]]['seq'])} bases)",
        key="studio_read_choose_v3",
    )
    return reads[chosen[0]]["seq"]


def render_teaching_studio(report, reads=None):
    st.markdown("### 🎬 Algorithm Teaching Studio — Moving DNA Search")
    st.caption("Compare → Shift/Fallback → Detect complete match. Pattern moves under the DNA. "
               "Play/Pause, step controls, a growing found list and exact-match notifications all run inside the visualization.")
    methods = report["methods"]
    targets = [Target(**x) for x in report["targets"]]
    c1, c2 = st.columns(2)
    with c1:
        method = st.selectbox("Algorithm to animate", methods, key="studio_method_v3")
    with c2:
        tgt_idx = st.selectbox("Pattern to trace (single-pattern methods)", list(range(len(targets))),
                               format_func=lambda i: f"{targets[i].name} ({targets[i].strand}) · {targets[i].sequence}",
                               key="studio_target_v3")
    effective = Matcher(method, targets).delegate.method if method == "Adaptive (ours)" else method
    trace_targets = targets if effective == "Aho-Corasick" else [targets[tgt_idx]]
    st.info(EXPLAIN[method] + (f"  Selected underlying algorithm: **{effective}**." if method == "Adaptive (ours)" else ""))

    input_mode = st.radio("Choose DNA for the animation", ["Real saved benchmark read", "Editable DNA example"],
                          horizontal=True, key="studio_input_mode_v3")
    default = report.get("example_read", "")
    if input_mode == "Real saved benchmark read":
        selected = _record_selector(reads, trace_targets, effective)
        if selected is not None:
            default = selected
        elif reads:
            st.caption("No matching read found in the loaded sample. The first read is shown instead.")
            default = reads[0]["seq"]
        st.code(f"Read length: {len(default)} bp | First 80 bases: {default[:80]}", language="text")
        seq = default
    else:
        examples = {
            "Overlapping — ACACACAC": "ACACACAC",
            "Repetitive / KMP fallback — AAAAA...C": "A" * 35 + "C",
            "Two exact hits — CACCCGCACCCG": "CACCCGCACCCG",
            "Use my own DNA": default,
        }
        preset = st.selectbox("Teaching example", list(examples), key="studio_example_v3")
        seq = st.text_area("DNA sequence for animation", value=examples[preset], height=100,
                           key=f"studio_dna_{preset.replace(' ','_')}").strip().upper()
    if not seq or set(seq) - set("ACGTRYSWKMBDHVN"):
        st.error("Enter DNA containing valid IUPAC nucleotide letters.")
        return
    if len(seq) > 10000:
        st.warning("This teaching player is intended for one manageable read (up to 10,000 bases). "
                   "The dataset benchmark still scans all selected reads. Choose a shorter read or select a region for animation.")
        return
    cap = min(100_000, max(5000, len(seq) * (max(len(t.sequence) for t in trace_targets) + 3) + 100))
    trace = build_teaching_trace(seq, trace_targets, effective, max_events=cap)
    if method == "Adaptive (ours)":
        trace["algorithm"] = method
        trace["strategy"] = f"Selected {effective}"
    if not trace["complete"]:
        st.warning("The teaching trace reached the 100,000-event browser safeguard. "
                   "This does not truncate your full dataset benchmark. Reduce sequence length to animate completely.")
    a, b, c = st.columns(3)
    a.metric("DNA bases in player", f"{len(seq):,}")
    b.metric("Playback events", f"{len(trace['events']):,}")
    c.metric("Exact hits in this read", str(len(trace["hits"])))
    components.html(
        _trace_player_html({"dna": seq, "trace": trace,
                            "targets": [{"name": t.name, "strand": t.strand, "sequence": t.sequence} for t in trace_targets],
                            "explain": EXPLAIN.get(effective, "")}),
        height=820, scrolling=False,
    )
    if effective == "Aho-Corasick":
        with st.expander("🔍 View actual Aho–Corasick trie and failure links"):
            nodes = Matcher("Aho-Corasick", trace_targets).trie
            lines = ['digraph Aho {rankdir=LR;node [shape=circle];']
            cap_nodes = min(len(nodes), 65)
            for idx in range(cap_nodes):
                lines.append(f'{idx} [shape={"doublecircle" if nodes[idx]["out"] else "circle"}];')
                for char, child in nodes[idx]["next"].items():
                    if child < cap_nodes:
                        lines.append(f'{idx} -> {child} [label="{char}"];')
                if idx and nodes[idx]["fail"] and nodes[idx]["fail"] < cap_nodes:
                    lines.append(f'{idx} -> {nodes[idx]["fail"]} [style=dashed,color=red,label="fail"];')
            lines.append('}')
            st.graphviz_chart(''.join(lines), use_container_width=True)
            st.caption(f"Showing {cap_nodes}/{len(nodes)} nodes; dashed red arrows are failure links.")
    if effective == "KMP":
        with st.expander("📐 Explain the KMP LPS table"):
            pattern = trace_targets[0].sequence
            st.dataframe([{"index": i, "nucleotide": b, "LPS": n} for i, (b, n) in enumerate(zip(pattern, lps_table(pattern)))], hide_index=True)
            st.caption("LPS[j−1] determines which prefix to reuse after a mismatch. Watch for FALLBACK frames in the player.")
    st.caption("Scientific scope: exact read-level motif positions, not chromosome coordinates or confirmed cancer variants. "
               "Python str.find cannot expose its native comparisons. Playback is one selected read; benchmark is all sampled reads.")
