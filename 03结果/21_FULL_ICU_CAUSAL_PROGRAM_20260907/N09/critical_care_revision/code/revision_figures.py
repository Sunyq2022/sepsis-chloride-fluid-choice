"""Supplementary diagnostic figure from aggregate revision outputs only."""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

REV=Path(__file__).resolve().parents[1]
d=pd.read_csv(REV/'results/propensity_distribution.csv')
plt.rcParams.update({'font.family':'Times New Roman','font.size':11,'text.color':'black',
 'axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.spines.top':False,
 'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,2,figsize=(8.0,3.6),sharex=True,sharey=True)
for ax,col,title in zip(axs,['proportion','weighted_proportion'],['A  Before weighting','B  After overlap weighting']):
    for arm,color,style in [('NS','#2678B2','-'),('LR','#D84A55','--')]:
        x=d[d.arm.eq(arm)]
        edges=list(x.left)+[float(x.right.iloc[-1])]
        ax.stairs(x[col],edges,color=color,linestyle=style,linewidth=1.8,label='Saline' if arm=='NS' else 'Lactated Ringer')
    ax.set(xlim=(0,1),ylim=(0,.26),xlabel='Propensity score for lactated Ringer')
    ax.set_title(title,loc='left',fontweight='bold',fontsize=12,pad=12)
    ax.yaxis.set_major_formatter(PercentFormatter(1,decimals=0));ax.legend(frameon=False,fontsize=9,loc='upper right')
axs[0].set_ylabel('Proportion within treatment group')
fig.subplots_adjust(left=.09,right=.985,bottom=.20,top=.86,wspace=.15)
for ext in ['png','pdf','svg']:
    fig.savefig(REV/'figures'/('Figure_S1_propensity_distributions.'+ext),dpi=300,facecolor='white')
plt.close(fig)
