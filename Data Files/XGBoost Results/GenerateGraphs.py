import pandas as pd
import matplotlib.pyplot as plt

font = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 9}
font1 = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 11}

font2 = {'family' : 'serif',
        'weight' : 'bold',
        'size'   : 13}

plt.rc('font', **font)
plt.rcParams['xtick.major.size'] = 3
plt.rcParams['xtick.major.width'] = 3
plt.rcParams['ytick.major.size'] = 3
plt.rcParams['ytick.major.width'] = 3
# plt.rcParams['xtick.minor.size'] = 10
# plt.rcParams['xtick.minor.width'] = 2

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, sharey=True, figsize=(13, 5))

###  LEARNING RATE = 0.01
input_file1 = 'CoulombMatrix.xlsx'
input_file2 = 'Mordred.xlsx'
input_file3 = 'Morgan.xlsx'
input_file4 = 'MACCS.xlsx'
input_file5 = 'Descriptor.xlsx'


CM = pd.read_excel(input_file1)
Mordred = pd.read_excel(input_file2)
Morgan = pd.read_excel(input_file3)
MACCS = pd.read_excel(input_file4)
Descriptor = pd.read_excel(input_file5)


exp_boiling = CM['Observed']
pred_boiling = CM['Predicted']

exp_boiling_wr_ind = CM['Difference'] > 30


ax1.plot([exp_boiling.min(), exp_boiling.max()], [exp_boiling.min(), exp_boiling.max()], 'k--', lw = 1)
ax1.scatter(exp_boiling[~exp_boiling_wr_ind], pred_boiling[~exp_boiling_wr_ind], color='b', marker='o')
ax1.scatter(exp_boiling[exp_boiling_wr_ind], pred_boiling[exp_boiling_wr_ind], color='r', marker='o')

# Hide the right and top spines
ax1.spines[['right', 'top']].set_visible(False)
ax1.spines[['left', 'bottom']].set_linewidth(2)

ax1.set_title('Coulomb Matrix', fontdict=font2)



exp_boiling = Mordred['Observed']
pred_boiling = Mordred['Predicted']

exp_boiling_wr_ind = Mordred['Difference'] > 30


ax2.plot([exp_boiling.min(), exp_boiling.max()], [exp_boiling.min(), exp_boiling.max()], 'k--', lw = 1)
ax2.scatter(exp_boiling[~exp_boiling_wr_ind], pred_boiling[~exp_boiling_wr_ind], color='b', marker='o')
ax2.scatter(exp_boiling[exp_boiling_wr_ind], pred_boiling[exp_boiling_wr_ind], color='r', marker='o')

# Hide the right and top spines
ax2.spines[['right', 'top']].set_visible(False)
ax2.spines[['left', 'bottom']].set_linewidth(2)

ax2.set_title('Mordred', fontdict=font2)






exp_boiling = Morgan['Observed']
pred_boiling = Morgan['Predicted']

exp_boiling_wr_ind = Morgan['Difference'] > 30


ax3.plot([exp_boiling.min(), exp_boiling.max()], [exp_boiling.min(), exp_boiling.max()], 'k--', lw = 1)
ax3.scatter(exp_boiling[~exp_boiling_wr_ind], pred_boiling[~exp_boiling_wr_ind], color='b', marker='o')
ax3.scatter(exp_boiling[exp_boiling_wr_ind], pred_boiling[exp_boiling_wr_ind], color='r', marker='o')

# Hide the right and top spines
ax3.spines[['right', 'top']].set_visible(False)
ax3.spines[['left', 'bottom']].set_linewidth(2)

ax3.set_title('Morgan', fontdict=font2)






exp_boiling = MACCS['Observed']
pred_boiling = MACCS['Predicted']

exp_boiling_wr_ind = MACCS['Difference'] > 30


ax4.plot([exp_boiling.min(), exp_boiling.max()], [exp_boiling.min(), exp_boiling.max()], 'k--', lw = 1)
ax4.scatter(exp_boiling[~exp_boiling_wr_ind], pred_boiling[~exp_boiling_wr_ind], color='b', marker='o')
ax4.scatter(exp_boiling[exp_boiling_wr_ind], pred_boiling[exp_boiling_wr_ind], color='r', marker='o')

# Hide the right and top spines
ax4.spines[['right', 'top']].set_visible(False)
ax4.spines[['left', 'bottom']].set_linewidth(2)

ax4.set_title('MACCS', fontdict=font2)




exp_boiling = Descriptor['Observed']
pred_boiling = Descriptor['Predicted']

exp_boiling_wr_ind = Descriptor['Difference'] > 30


ax5.plot([exp_boiling.min(), exp_boiling.max()], [exp_boiling.min(), exp_boiling.max()], 'k--', lw = 1)
ax5.scatter(exp_boiling[~exp_boiling_wr_ind], pred_boiling[~exp_boiling_wr_ind], color='b', marker='o')
ax5.scatter(exp_boiling[exp_boiling_wr_ind], pred_boiling[exp_boiling_wr_ind], color='r', marker='o')

# Hide the right and top spines
ax5.spines[['right', 'top']].set_visible(False)
ax5.spines[['left', 'bottom']].set_linewidth(2)

ax5.set_title('Descriptor', fontdict=font2)

fig.text(0.5, 0.04, 'Observed Melting Points $^{\circ}$C', ha='center')
fig.text(0.08, 0.5, 'Predicted Melting Points $^{\circ}$C', va='center', rotation='vertical')

plt.subplots_adjust(wspace=0, hspace=0)
plt.show()
