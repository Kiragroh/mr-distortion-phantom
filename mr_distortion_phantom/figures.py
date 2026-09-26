"""Rebuild manuscript plots and editable tables from the released measurements."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Circle
from .core import evaluate, metrics, sphere_coverage

BLUE='#008ac9';PURPLE='#8064a2';INK='#243357';ORANGE='#c66b27'
DATA=Path(__file__).parent/'data'

def load():
    groups={}
    with (DATA/'landmarks.csv').open() as f:
        for row in csv.DictReader(f):groups.setdefault((row['protocol'],row['state']),[]).append(row)
    output={}
    for key,rows in groups.items():
        ct=np.array([[float(r['ct_'+a+'_mm']) for a in 'LPS'] for r in rows]);mr=np.array([[float(r['mr_'+a+'_mm']) for a in 'LPS'] for r in rows]);fit=np.array([r['role']=='registration' for r in rows])
        stat,vec,_,_=evaluate(ct,mr,fit)
        output[key]=dict(stats=stat,errors=np.linalg.norm(vec[~fit],axis=1),fit=metrics(np.linalg.norm(vec[fit],axis=1)),ct=ct,mask=fit,rows=rows)
    return output

def save(fig,path):
    for ext in ['svg','pdf','png']:fig.savefig(path.with_suffix('.'+ext),dpi=300,bbox_inches='tight',facecolor='white',metadata={'Creator':'mr-distortion-phantom'} if ext in ['pdf','svg'] else None)
    plt.close(fig)

def table(out,name,headers,rows,caption):
    with (out/(name+'.csv')).open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f);w.writerow(headers);w.writerows(rows)
    md='**'+name.replace('_',' ')+'.** '+caption+'\n\n| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'
    md+='\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)+'\n'
    (out/(name+'.md')).write_text(md,encoding='utf8')
    import html
    ht='<table><caption>'+html.escape(caption)+'</caption><thead><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in headers)+'</tr></thead><tbody>'
    ht+=''.join('<tr>'+''.join('<td>'+html.escape(str(x))+'</td>' for x in r)+'</tr>' for r in rows)+'</tbody></table>'
    (out/(name+'.html')).write_text(ht,encoding='utf8')
    def tex(x):return str(x).replace('_',r'\_').replace('%',r'\%').replace('Δ',r'$\Delta$').replace('−','-').replace('≤',r'$\leq$')
    latex='\\begin{table}\n\\caption{'+tex(caption)+'}\n\\centering\n\\begin{tabular}{'+'l'*len(headers)+'}\n\\hline\n'+' & '.join(map(tex,headers))+r' \\'+'\n\\hline\n'
    latex+='\n'.join(' & '.join(map(tex,r))+r' \\' for r in rows)+'\n\\hline\n\\end{tabular}\n\\end{table}\n'
    (out/(name+'.tex')).write_text(latex,encoding='utf8')

def boxes(ax,arrays,positions,color,width=.25):
    bp=ax.boxplot(arrays,positions=positions,vert=False,widths=width,whis=(5,95),showfliers=False,patch_artist=True,manage_ticks=False)
    for b in bp['boxes']:b.set(facecolor=color,alpha=.25,edgecolor=color)
    for k in ['whiskers','caps','medians']:
        for b in bp[k]:b.set(color=color,linewidth=1.6)

def build(out):
    out=Path(out)
    for d in ['figures','tables','supplement','docs']:(out/d).mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'text.color':INK,'axes.labelcolor':INK,'svg.fonttype':'none','pdf.fonttype':42})
    g=load();meta=json.loads((DATA/'protocols.json').read_text());rep=json.loads((DATA/'report_distributions.json').read_text());val=json.loads((DATA/'technical_validation.json').read_text())
    pairs=[p['protocol'] for p in meta if p['paired']];codes=[p['protocol'] for p in meta];captions={}
    # Figure 1: measured CT coordinates and disjoint registration/evaluation roles.
    fig=plt.figure(figsize=(11,5.2),layout='constrained');gs=fig.add_gridspec(1,2,width_ratios=[1.1,1]);ax=fig.add_subplot(gs[0]);ax.axis('off')
    steps=[('1  Native phantom images','Planning CT + original/corrected MRI'),('2  Detect and index grid nodes','DoG peaks + subvoxel localization; no report input'),('3  Rigid alignment, separately per MRI','91 central nodes; nominal radius ≤30 mm'),('4  Evaluate identical outer nodes','757 held-out nodes; paired protocol distributions')]
    for i,(h,t) in enumerate(steps):
        y=.91-i*.245;ax.text(.03,y,h,fontsize=13,weight='bold',va='top',color=BLUE);ax.text(.03,y-.075,t,fontsize=10,va='top')
        if i<3:ax.annotate('',(.1,y-.205),(.1,y-.13),arrowprops={'arrowstyle':'->','color':PURPLE})
    ax=fig.add_subplot(gs[1],projection='3d');v=g[('P02','original')];xyz=v['ct'];mask=v['mask'];ax.scatter(*xyz[~mask].T,c=BLUE,s=6,alpha=.3,label='Evaluation (n=757)');ax.scatter(*xyz[mask].T,c=PURPLE,s=13,label='Registration (n=91)');ax.set(xlabel='Centred CT L [mm]',ylabel='Centred CT P [mm]',zlabel='Centred CT S [mm]');ax.legend(loc='upper left',fontsize=8);ax.view_init(23,32)
    save(fig,out/'figures/Figure_1_Workflow')
    captions['Figure 1']='CT-referenced workflow and measurement geometry. Grid coordinates are independently detected from native phantom images; report values do not enter the analysis. Each MRI is rigidly aligned to the measured CT using 91 central nodes, selected by a nominal lattice radius of 30 mm. Residual displacement is evaluated at the remaining 757 common nodes. The 3D view shows measured, centred CT coordinates for P02. The nominal lattice assigns correspondence and radius only; it is not a replacement for measured CT or the manufacturer CAD reference.'
    # Figure 2: distributions and paired P95; identical order and n.
    fig,(ax,b)=plt.subplots(1,2,figsize=(11,5.7),gridspec_kw={'width_ratios':[2.6,1]},layout='constrained');ys=np.arange(len(pairs))
    for state,off,col in [('original',-.17,PURPLE),('corrected',.17,BLUE)]:boxes(ax,[g[(p,state)]['errors'] for p in pairs],ys+off,col)
    ax.set(yticks=ys,yticklabels=pairs,xlabel='CT-referenced residual displacement [mm]',title='A  Per-protocol distributions');ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    for i,p in enumerate(pairs):
        a=g[(p,'original')]['stats']['p95_mm'];c=g[(p,'corrected')]['stats']['p95_mm'];b.plot([a,c],[i,i],color='#a3aaba');b.scatter(a,i,color=PURPLE,s=35);b.scatter(c,i,color=BLUE,s=35);b.text(max(a,c)+.025,i,f'{c-a:+.3f}',va='center',fontsize=9)
    b.set(yticks=ys,yticklabels=pairs,xlabel='P95 [mm]; label = ΔP95',title='B  Paired P95');b.set_xlim(0,max(x['stats']['p95_mm'] for x in g.values())+.22);b.invert_yaxis();b.grid(axis='x',alpha=.15)
    fig.legend(handles=[Patch(color=PURPLE,label='Original MRI'),Patch(color=BLUE,label='Corrected MRI')],loc='outside lower center',ncol=2,frameon=False)
    save(fig,out/'figures/Figure_2_Paired_distributions')
    captions['Figure 2']='Protocol-specific CT-referenced residuals for seven original–corrected MRI pairs. Each box contains the same 757 held-out nodes within its protocol: median, interquartile range and 5th–95th percentile whiskers (points outside whiskers are not displayed). Panel B shows P95 before and after correction; labels give ΔP95 = P95(corrected) − P95(original), in mm. Quantiles use NumPy linear interpolation. Separate rigid fits remove central rigid misalignment from each series. Nodes are spatially dependent and are not independent acquisition replicates; no pooled significance test is implied.'
    # Figure 3: aggregate differences, all 16 series.
    keys=sorted(g);fig,(ax,b)=plt.subplots(1,2,figsize=(11,6.5),sharey=True,layout='constrained');ys=np.arange(len(keys));ctrep={(r['protocol'],r['state']):r for r in rep if r['reference']=='CT'}
    for axis,metric,title in [(ax,'mean_mm','A  Mean residual'),(b,'p95_mm','B  95th percentile')]:
        for i,k in enumerate(keys):
            r=ctrep[k];rv=r['mean_mm'] if metric=='mean_mm' else r['quantiles_mm']['P95'];own=g[k]['stats'][metric]
            axis.plot([rv,own],[i,i],color='#c7cbd3');axis.scatter(rv,i,c=ORANGE,marker='s',s=27);axis.scatter(own,i,c=BLUE,s=30)
        axis.set(xlabel='Residual [mm]',title=title);axis.grid(axis='x',alpha=.15)
    ax.set(yticks=ys,yticklabels=[p+' '+('C' if s=='corrected' else 'O') for p,s in keys]);ax.invert_yaxis();fig.legend(handles=[Patch(color=BLUE,label='Independent calculation (757 nodes)'),Patch(color=ORANGE,label='CT-reference report (858 nodes)')],loc='outside lower center',ncol=2,frameon=False)
    save(fig,out/'figures/Figure_3_Report_comparison')
    captions['Figure 3']='Descriptive comparison of aggregate residual metrics from the independent calculation and CT-referenced manufacturer reports across all 16 eligible MRI series (O: original; C: corrected). Series and row order are identical in both panels. Independent calculations use 757 held-out nodes after central rigid alignment; reports contain 858 nodes with unmatched identities and a different analysis convention. Therefore, connecting lines show aggregate differences, not nodewise agreement or equivalence. Report values were recovered from vector graphics and checked against printed report summaries. The two largest mean differences occur in the P08 and P09 original series.'
    # Figure 4: technical tests, explicitly show outlier tails as well as P95.
    fig,(ax,b)=plt.subplots(1,2,figsize=(11,5.3),gridspec_kw={'width_ratios':[1.4,1]},layout='constrained');ts=sorted(val['translations'],key=lambda x:x['protocol']);ys=np.arange(len(ts))
    ax.scatter([t['p95_mm'] for t in ts],ys,color=BLUE,label='P95',s=40);ax.scatter([t['max_mm'] for t in ts],ys,color=ORANGE,label='Maximum',marker='x',s=45);ax.set(yticks=ys,yticklabels=[t['protocol'] for t in ts],xscale='log',xlabel='Translation recovery error [mm] (log scale)',title='A  Known voxel translations');ax.invert_yaxis();ax.grid(axis='x',alpha=.15);ax.legend(frameon=False)
    sig=[x['sigma_mm'] for x in val['sensitivity']]+[1.4];bef=[x['before_mean_mm'] for x in val['sensitivity']]+[g[('P08','original')]['stats']['mean_mm']];aft=[x['after_mean_mm'] for x in val['sensitivity']]+[g[('P08','corrected')]['stats']['mean_mm']];order=np.argsort(sig)
    for z,c,l in [(bef,PURPLE,'Original'),(aft,BLUE,'Corrected')]:b.plot(np.array(sig)[order],np.array(z)[order],'o-',color=c,label=l)
    b.set(xlabel='Detection σ [mm]',ylabel='Mean CT residual [mm]',title='B  P08 detector sensitivity',ylim=(0,.4));b.grid(alpha=.15);b.legend(frameon=False)
    save(fig,out/'figures/Figure_4_Technical_verification')
    captions['Figure 4']='Technical verification and parameter sensitivity. A: Original MRI volumes were translated by (+0.35, −0.45, +0.25) mm in LPS using cubic interpolation, followed by landmark redetection. P95 recovery error was below 0.094 mm in all nine protocols; maxima reached 1.789 mm, illustrating the need to inspect rare localization/correspondence outliers. These tests assess numerical equivariance, not absolute physical accuracy. B: For P08, the mean corrected residual remained below the original residual at detection σ values of 1.2, 1.4 and 1.6 mm (second DoG scale fixed at 3.5 mm). Sensitivity analyses do not resolve the remaining report discrepancy.'
    # Supplement: CT versus CAD report distributions; separate protocol rows, no pooling.
    fig,axes=plt.subplots(1,2,figsize=(11,6),sharex=True,sharey=True,layout='constrained');ys=np.arange(9)
    for ax,state,title in zip(axes,['original','corrected'],['Original MRI','Corrected MRI']):
        for ref,off,col in [('CT',-.17,BLUE),('CAD',.17,ORANGE)]:
            selected=[r for r in rep if r['state']==state and r['reference']==ref];mapping={r['protocol']:r for r in selected};present=[(i,p) for i,p in enumerate(codes) if p in mapping];boxes(ax,[mapping[p]['errors_mm'] for i,p in present],np.array([i for i,p in present])+off,col)
        for i,p in enumerate(codes):
            if state=='corrected' and p not in pairs:ax.text(.1,i,'No eligible corrected acquisition',fontsize=8,color='#777')
        ax.set(title=title,xlabel='Report residual [mm]',yticks=ys,yticklabels=codes);ax.grid(axis='x',alpha=.15)
    axes[0].invert_yaxis();fig.legend(handles=[Patch(color=BLUE,label='Planning CT (858 nodes)'),Patch(color=ORANGE,label='CAD (859 nodes)')],loc='outside lower center',ncol=2,frameon=False);save(fig,out/'supplement/Figure_S1_CT_CAD_reports')
    captions['Figure S1']='Report-derived residual distributions for each protocol against the planning CT and manufacturer CAD references. Boxes show median and interquartile range, with 5th–95th percentile whiskers; CT reports contain 858 nodes, CAD reports 859. Original and corrected panels retain identical protocol order; missing corrected acquisitions remain explicitly visible. These are report-derived distributions, not independent CAD measurements. Similarity of their summaries does not demonstrate interchangeability, which would require matched coordinates and a prespecified equivalence margin.'
    # Synthetic coverage only, deliberately independent of patient examples.
    fig,(ax,b)=plt.subplots(1,2,figsize=(11,4.8),gridspec_kw={'width_ratios':[1.6,1]},layout='constrained');D=np.linspace(3,40,400)
    for d,c in [(.5,BLUE),(1,PURPLE),(2,ORANGE)]:ax.plot(D,100*sphere_coverage(D,d),label=f'{d:g} mm shift',color=c,lw=2)
    ax.set(xlabel='Sphere diameter [mm]',ylabel='Geometric overlap [%]',ylim=(0,101),title='Synthetic equal-sphere model');ax.grid(alpha=.15);ax.legend(frameon=False)
    b.add_patch(Circle((0,0),5,color=BLUE,alpha=.25));b.add_patch(Circle((2,0),5,fill=False,edgecolor=ORANGE,lw=2));b.annotate('',(2,0),(0,0),arrowprops={'arrowstyle':'<->','color':INK});b.text(1,.6,'2 mm',ha='center');b.set(xlim=(-6,8),ylim=(-6,6),aspect='equal',title='10 mm diameter; 70.4% overlap');b.axis('off');save(fig,out/'supplement/Figure_S2_Synthetic_coverage')
    captions['Figure S2']='Synthetic illustration of size-dependent geometric overlap. Two spheres with identical radius R and centre separation d have overlap fraction 1 − 3d/(4R) + d³/(16R³) for 0 ≤ d ≤ 2R, and zero beyond. Curves show chosen displacements of 0.5, 1 and 2 mm; the diagram is a 2D section of the 3D model. This idealized geometric quantity is not measured target coverage, a dose-volume metric, or evidence of clinical benefit. No patient structures or observed patient displacements were used.'
    # Editable tables.
    t1=[]
    for p in meta:t1.append([p['protocol'],p['scanner'],p['manufacturer'],f"{p['field_T']:g}",f"{p['TR_ms']:.2f}",f"{p['TE_ms']:.2f}",f"{p['bandwidth_Hz_px']:.1f}",' × '.join(f'{v:.3f}' for v in p['voxel_mm']),'Yes' if p['paired'] else 'No'])
    table(out/'tables','Table_1_Acquisitions',['Protocol','Scanner','Manufacturer','B0 [T]','TR [ms]','TE [ms]','BW [Hz/pixel]','Voxel z/y/x [mm]','Paired'],t1,'Acquisition inventory: nine original protocols from five scanner groups, with seven eligible corrected counterparts (16 MRI series total). One planning CT is the shared reference. Scanner codes do not denote independent centres; DICOM TR definitions may differ between sequence families and manufacturers.')
    t2=[]
    for p in pairs:
        a=g[(p,'original')]['stats'];c=g[(p,'corrected')]['stats'];t2.append([p,757,f"{a['median_mm']:.3f}",f"{c['median_mm']:.3f}",f"{a['p95_mm']:.3f}",f"{c['p95_mm']:.3f}",f"{c['p95_mm']-a['p95_mm']:+.3f}"])
    table(out/'tables','Table_2_Paired_results',['Protocol','n','Median O [mm]','Median C [mm]','P95 O [mm]','P95 C [mm]','ΔP95 [mm]'],t2,'Paired CT-referenced residuals at 757 identical evaluation nodes per protocol after separate central rigid fits. O: original; C: corrected. ΔP95 is the difference between two distribution quantiles, not the 95th percentile of nodewise changes. Endpoints were selected post hoc and are descriptive.')
    t3=[['Public coordinate reproduction','16 series × 848 nodes','All summary metrics reproduced within 1e-9 mm','Checks numeric reproducibility, not new accuracy evidence'],['Known voxel translations','9 original MRI volumes',f"P95 {min(t['p95_mm'] for t in ts):.3f}–{max(t['p95_mm'] for t in ts):.3f} mm; maximum {max(t['max_mm'] for t in ts):.3f} mm",'Interpolation / localization / correspondence test'],['Detector sensitivity','P08; σ 1.2–1.6 mm','Mean reduction 0.119–0.128 mm','One-protocol sensitivity only'],['Independent vs report mean','16 unmatched-node comparisons','14 absolute differences <0.036 mm; two ≈0.105 mm','Descriptive comparison; no equivalence claim']]
    table(out/'tables','Table_3_Verification',['Test','Scope','Result','Interpretation'],t3,'Verification layers and their scope. Repeated/repositioned acquisitions and independent physical ground truth are not available in this dataset.')
    allrows=[]
    for k in keys:
        a=g[k]['stats'];r=ctrep[k];allrows.append([*k,a['n'],*[f'{a[x]:.6f}' for x in ['mean_mm','median_mm','p95_mm','max_mm']],f"{g[k]['fit']['mean_mm']:.6f}",r['n'],f"{r['mean_mm']:.6f}",f"{r['quantiles_mm']['P95']:.6f}",f"{a['mean_mm']-r['mean_mm']:+.6f}"])
    table(out/'supplement','Table_S1_All_series',['Protocol','State','n eval','Own mean','Own median','Own P95','Own max','Central fit mean','n report','Report mean','Report P95','Mean difference'],allrows,'All-series results in mm, including central-fit residuals and aggregate CT-report comparisons. Node identities differ between methods; report quantiles originate from extracted vector points.')
    table(out/'supplement','Table_S2_Coverage',['Diameter [mm]','Shift [mm]','Geometric overlap [%]'],[[D,d,f'{100*sphere_coverage(D,d):.2f}'] for D in [5,10,20,30] for d in [.5,1,2]],captions['Figure S2'])
    (out/'docs/CAPTIONS.md').write_text('# Figure captions\n\n'+'\n\n'.join('## '+k+'\n\n'+v for k,v in captions.items())+'\n\nPlots were created using conventional numerical analysis and Matplotlib with AI-assisted code development. No generative image model was used. The authors must verify and disclose assistance according to the journal policy.\n',encoding='utf8')
    (out/'docs/captions.json').write_text(json.dumps(captions,ensure_ascii=False,indent=2),encoding='utf8')
    summaries={'pairs':len(pairs),'delta_p95_mm':{p:g[(p,'corrected')]['stats']['p95_mm']-g[(p,'original')]['stats']['p95_mm'] for p in pairs},'all_mean_difference_mm':{p+'_'+s:g[(p,s)]['stats']['mean_mm']-ctrep[(p,s)]['mean_mm'] for p,s in keys},'translation_p95_max_mm':max(t['p95_mm'] for t in ts),'translation_max_mm':max(t['max_mm'] for t in ts)}
    (out/'docs/result_summary.json').write_text(json.dumps(summaries,indent=2))
    print('Created four main figures, two supplement figures, three main tables and two supplement tables.')

if __name__=='__main__':
    import sys
    build(Path(sys.argv[1]) if len(sys.argv)>1 else Path('.'))
