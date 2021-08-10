import pandas as pd

################################################################################
#
# This script generates a plot from the output of the modification_script.py.
# The plot is an interactive bokeh plot, and comes in the form of an .html file.
#
################################################################################

from bokeh.palettes import Dark2
from bokeh.plotting import figure, output_file, show
from bokeh.models import Span

df_a_vs_i = pd.read_csv('/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/19_6mer_comp_reprise/A_vs_I/results/out_nanocompore_shift_stats.tsv', sep='\t')
print(df_a_vs_i.shape)
df_c_vs_i = pd.read_csv('/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/19_6mer_comp_reprise/C_vs_I/results/out_nanocompore_shift_stats.tsv', sep='\t')
print(df_c_vs_i.shape)
df_g_vs_i = pd.read_csv('/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/19_6mer_comp_reprise/G_vs_I/results/out_nanocompore_shift_stats.tsv', sep='\t')
print(df_g_vs_i.shape)
df_t_vs_i = pd.read_csv('/export/valenfs/data/processed_data/MinION/inosine_project/20201016_max_DNA_Inosine/second_analysis/19_6mer_comp_reprise/T_vs_I/results/out_nanocompore_shift_stats.tsv', sep='\t')
print(df_t_vs_i.shape)

df_compiled = pd.DataFrame(index=range(0,64), columns=['A_intensity', 'C_intensity', 'G_intensity', 'T_intensity', 'I_A_ref_intensity', 'I_C_ref_intensity','I_G_ref_intensity','I_T_ref_intensity'], dtype='float')
print(df_compiled.shape)
df_compiled['A_intensity'] = df_a_vs_i['c1_mean_intensity']
df_compiled['C_intensity'] = df_c_vs_i['c1_mean_intensity']
df_compiled['G_intensity'] = df_g_vs_i['c1_mean_intensity']
df_compiled['T_intensity'] = df_t_vs_i['c1_mean_intensity']
df_compiled['I_A_ref_intensity'] = df_a_vs_i['c2_mean_intensity']
df_compiled['I_C_ref_intensity'] = df_c_vs_i['c2_mean_intensity']
df_compiled['I_G_ref_intensity'] = df_g_vs_i['c2_mean_intensity']
df_compiled['I_T_ref_intensity'] = df_t_vs_i['c2_mean_intensity']

print(df_compiled.head())

p = figure(plot_width=800, plot_height=250)
p.title.text = 'Click on legend entries to hide the corresponding lines'
p.ygrid.grid_line_color = None

vline1 = Span(location=15, dimension="height", line_color='black', line_width=0.5)
vline2 = Span(location=29, dimension="height", line_color='black', line_width=0.5)
vline3 = Span(location=44, dimension="height", line_color='black', line_width=0.5)

for col, color in zip(df_compiled.columns,Dark2[8]):
	print(color)
	p.line(df_compiled.index, df_compiled[col], line_width=1.2, color=color, alpha=0.8,legend_label=col)
p.legend.location = "top_left"
p.legend.click_policy = "hide"

output_file("8_curve_plot.html", title="interactive_legend.py example")
p.renderers.extend([vline1, vline2, vline3])
show(p)
