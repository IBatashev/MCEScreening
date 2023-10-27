import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('Tc_data.csv', index_col=0, sep=',')

x = df['magF_u'].to_numpy()
y = df['Tc'].to_numpy()
names = df['pretty_formula'].to_numpy()

# norm = plt.Normalize(1,4)
# cmap = plt.cm.RdYlGn

fig, ax = plt.subplots(figsize=(10, 10))
sc = plt.scatter(x, y, s=50, color='black')

annot = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points", fontsize=17,
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
annot.set_visible(False)


def update_annot(ind):
    pos = sc.get_offsets()[ind["ind"][0]]
    annot.xy = pos
    text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_alpha(0.4)


def hover(event):
    vis = annot.get_visible()
    if event.inaxes == ax:
        cont, ind = sc.contains(event)
        if cont:
            update_annot(ind)
            annot.set_visible(True)
            fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                fig.canvas.draw_idle()

fig.canvas.mpl_connect("motion_notify_event", hover)


plt.xlim(0.45, 2)
plt.ylim(0, 1000)

ax.tick_params(axis='both', which='major', labelsize=17)
plt.xlabel("Internal Field, T", fontsize=17)
plt.ylabel("$T_{C}$, K", fontsize=17)

plt.show()