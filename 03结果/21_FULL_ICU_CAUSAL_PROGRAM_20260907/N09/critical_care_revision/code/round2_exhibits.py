"""Current aggregate baseline summaries and a restrained forest display."""
import pandas as pd,numpy as np
import revision_analysis as a
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# Reuse the established response figure palette and typography without the docx dependency.
plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
R=a.REV;O=R/'results/round2'
d=pd.read_parquet(R/'work/private/round2/cohort.parquet')
dd,prep,x,ps,w=a.n.fit(d,a.m,False)
summaries=[]
for col in a.n.NUM:
    row={'variable':col};missing=[]
    for arm,label in [(0,'NS'),(1,'LR')]:
        z=d.loc[d.exposure.eq(arm),col];row[label]=f'{z.mean():.2f} ({z.std(ddof=0):.2f})';missing.append(str(z.isna().sum()))
    row['missing']=' / '.join(missing);summaries.append(row)
pd.DataFrame(summaries).to_csv(O/'covariate_summary.csv',index=False)
cat=[]
for col in a.n.CAT:
    for arm,label in [(0,'NS'),(1,'LR')]:
        for level,count in d.loc[d.exposure.eq(arm),col].value_counts(dropna=False).items():cat.append({'variable':col,'arm':label,'level':level,'n':count})
pd.DataFrame(cat).to_csv(O/'categorical_summary.csv',index=False)
hist=[];support=[]
for arm,label in [(0,'NS'),(1,'LR')]:
    mask=dd.exposure.to_numpy(int)==arm
    raw,edges=np.histogram(ps[mask],bins=np.linspace(0,1,21))
    weighted,_=np.histogram(ps[mask],bins=edges,weights=w[mask])
    for j in range(20):hist.append({'arm':label,'left':edges[j],'raw':raw[j]/mask.sum(),'weighted':weighted[j]/w[mask].sum()})
    support.append({'arm':label,'ps_min':ps[mask].min(),'ps_max':ps[mask].max(),'weight_min':w[mask].min(),'weight_max':w[mask].max()})
pd.DataFrame(hist).to_csv(O/'propensity_distribution.csv',index=False)
pd.DataFrame(support).to_csv(O/'propensity_support.csv',index=False)
fig,axes=plt.subplots(1,2,figsize=(8,2.8),sharey=True)
for ax,key,title in zip(axes,['raw','weighted'],['Before weighting','After overlap weighting']):
    for arm,label,color,ls in [(0,'NS','#2678B2','-'),(1,'LR','#D84A55','--')]:
        h=[r for r in hist if r['arm']==label]
        ax.stairs([r[key] for r in h],np.linspace(0,1,21),color=color,ls=ls,lw=1.5,label=label)
    ax.set(title=title,xlabel='Propensity for LR',xlim=(0,1));ax.legend(frameon=False)
axes[0].set_ylabel('Within-arm proportion');fig.tight_layout()
for ext in ['png','pdf','svg']:fig.savefig(R/'figures'/('Figure_round2_propensity.'+ext),dpi=300,bbox_inches='tight')
plt.close(fig)
z=d.lr_ml.gt(0).to_numpy(float);aa=d.exposure.to_numpy(int);sd=np.sqrt((z[aa==0].var()+z[aa==1].var())/2)
pd.DataFrame([{'variable':'Any prior LR','smd_after':abs(np.average(z[aa==1],weights=w[aa==1])-np.average(z[aa==0],weights=w[aa==0]))/sd}]).to_csv(O/'prior_lr_balance.csv',index=False)
r=pd.read_csv(O/'estimates.csv');r=r[r.outcome.eq('primary')&~r.analysis.eq('Hospital follow-up')]
fig,ax=plt.subplots(figsize=(8.7,4.8));fig.subplots_adjust(left=.32,right=.98,bottom=.16,top=.94)
for y,(_,row) in zip(range(len(r)-1,-1,-1),r.iterrows()):
    value,lo,hi=[row[k] for k in ['rd_ns_lr_pp','rd_ns_lr_low_pp','rd_ns_lr_high_pp']]
    color='#168C8C' if row.analysis in ['Absolute only','Rise only','Both'] else '#2678B2'
    ax.errorbar(value,y,xerr=[[value-lo],[hi-value]],fmt='o',color=color,ms=6,capsize=3,lw=1.4)
    ax.text(26,y,f'{value:.2f} ({lo:.2f}, {hi:.2f})',va='center',fontsize=10)
ax.axvline(0,color='black',ls='--',lw=.8)
ax.set(yticks=list(range(len(r)-1,-1,-1)),yticklabels=[f'{x.analysis.replace(">=", "≥")} (n = {x.n})' for x in r.itertuples()],xlim=(-15,47),ylim=(-.5,len(r)-.3),xticks=[-10,0,10,20],xlabel='Risk difference, saline − LR (percentage points)')
ax.spines['left'].set_visible(False);ax.spines['bottom'].set_bounds(-15,25);ax.tick_params(axis='y',length=0,pad=8)
ax.text(26,len(r)-.35,'RD (95% CI)',fontsize=10,weight='bold')
for ext in ['png','pdf','svg']:fig.savefig(R/'figures'/('Figure_round2_estimates.'+ext),dpi=300,facecolor='white',bbox_inches='tight',pad_inches=.12)
plt.close(fig)
