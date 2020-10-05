from os.path import dirname, join

import numpy as np
import pandas as pd

from bokeh.io import curdoc
from bokeh.layouts import column, layout
from bokeh.models import ColumnDataSource, Div, Select, Slider, TextInput, DataTable, TableColumn, NumberFormatter, CustomJS, Button
from bokeh.plotting import figure
from bokeh.models import HoverTool
from bokeh.models.widgets import HTMLTemplateFormatter

datalist = pd.read_csv('D:/MCES/MP/step1_success_sieved_out.csv')
datalist.fillna(0, inplace=True)  # just replace missing values with zero

datalist["color"] = np.where(datalist["Mel_full"] >= 0, "orange", "grey")
# datalist["alpha"] = np.where(datalist["Mel_full"] > 0, 0.9, 0.25)

# with open(join(dirname(__file__), "razzies-clean.csv")) as f:
#     razzies = f.read().splitlines()
# movies.loc[movies.imdbID.isin(razzies), "color"] = "purple"
# movies.loc[movies.imdbID.isin(razzies), "alpha"] = 0.9

axis_map = {
    "Magneto Elastic": "Mel_full",
    "Mag. El. a": "Mel_a",
    "Mag. El. b": "Mel_b",
    "Mag. El. c": "Mel_c",
    "Internal Field Calculated": "magF_u",
    "Internal Field Database": "mag_field",
    "Volume": "volume_cell",
}


desc = Div(text=open(join(dirname(__file__), "Title.html")).read(), sizing_mode="stretch_width", max_width=1400)
text = Div(text=open(join(dirname(__file__), "Description.html")).read(), sizing_mode="stretch_width", max_width=320, margin=(40, 40, 40, 40))
# Create Input controls
# reviews = Slider(title="Magneto Elastic", value=80, start=10, end=300, step=10)
# min_year = Slider(title="Internal Field", start=1940, end=2014, value=1970, step=1)
# max_year = Slider(title="End Year released", start=1940, end=2014, value=2014, step=1)
# oscars = Slider(title="Minimum number of Oscar wins", start=0, end=4, value=0, step=1)
# boxoffice = Slider(title="Dollars at Box Office (millions)", start=0, end=800, value=0, step=1)
# genre = Select(title="Genre", value="All", options=open(join(dirname(__file__), 'symmetries.txt')).read().split())
# director = TextInput(title="Contains Element")
# cast = TextInput(title="Compound name contains")

symType = Select(title="Symmetry", value="Any", options=open(join(dirname(__file__), 'symmetries.txt')).read().split())

magEl = Slider(title="Magneto Elastic", value=0, start=0, end=datalist['Mel_full'].max(), step=0.5)
intField = Slider(title="Internal Field", start=0, end=datalist['magF_u'].max(), value=0, step=0.05)

compName = TextInput(title="Compound name contains")
elemIn = TextInput(title="Contains Element(s)")
elemOut = TextInput(title="Does not contain Element(s)")

x_axis = Select(title="X Axis", options=sorted(axis_map.keys()), value="Magneto Elastic")
y_axis = Select(title="Y Axis", options=sorted(axis_map.keys()), value="Internal Field Calculated")

button = Button(label="Save Table", button_type="success", height=50, width=150, sizing_mode="fixed")



# Create Column Data Source that will be used by the plot
source = ColumnDataSource(data=dict(x=[],
                                    y=[],
                                    color=[],
                                    id=[],
                                    compound=[],
                                    spacegroup=[],
                                    polymorphs=[],
                                    mag_sites=[],
                                    mag_type=[],
                                    intFieldDB=[],
                                    intField=[],
                                    magEl=[],
                                    magEl_a=[],
                                    magEl_b=[],
                                    magEl_c=[],
                                    volume=[],
                                    symmetry=[])
                          )
source_back = ColumnDataSource(data=dict(x=[],
                                         y=[],
                                         compound=[],
                                         intField=[],
                                         magEl=[],
                                         symmetry=[])
                               )

TOOLTIPS=[
    ("Compound", "@compound"),
    ("Symmetry", "@symmetry"),
    ("Mel", "@magEl"),
    ("Field", "@intField")
]

hover = HoverTool(names=['needshover'], tooltips=TOOLTIPS)

p = figure(plot_height=550, plot_width=550, title="", toolbar_location="right", tools="pan,wheel_zoom,box_zoom,lasso_select,tap,reset,save", active_scroll="wheel_zoom", tooltips=TOOLTIPS, sizing_mode="fixed")

p.circle(x="x", y="y", source=source_back, size=4, color="black", line_color=None, fill_alpha=0.5)
p.circle(x="x", y="y", source=source, size=5, color="color", line_color=None, fill_alpha=0.9, name='needshover')

p.add_tools(hover)
p.toolbar.active_inspect = [hover]
p.axis.axis_label_text_font_size = "20pt"
p.title.text_font_size = "20pt"

p.xaxis.major_label_text_font_size = "20pt"
p.yaxis.major_label_text_font_size = "20pt"

s2 = ColumnDataSource(data=dict(id=[], compound=[], symmetry=[]))

columns = [
    TableColumn(field="id", title="ID"),
    TableColumn(field="compound", title="Compound"),
    TableColumn(field="symmetry", title="Symmetry"),
    TableColumn(field="spacegroup", title="Space group"),
    TableColumn(field="polymorphs", title="Polymorphs"),
    TableColumn(field="mag_type", title="Ordering"),
    TableColumn(field="mag_sites", title="Sites"),
    TableColumn(field="intFieldDB", title="Internal Field DB", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="intField", title="Internal Field", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="volume", title="Volume", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="magEl", title="Magneto Elastic", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="magEl_a", title="Mag.El. a", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="magEl_b", title="Mag.El. b", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="magEl_c", title="Mag.El. c", formatter=NumberFormatter(format="0.00")),
]
data_table = DataTable(source=s2, columns=columns, width=600)

source.selected.js_on_change('indices', CustomJS(args=dict(s1=source, s2=s2), code="""
        var inds = cb_obj.indices;
        var d1 = s1.data;
        var d2 = s2.data;
        d2['id'] = []
        d2['compound'] = []
        d2['symmetry'] = []
        d2['spacegroup'] = []
        d2['polymorphs'] = []
        d2['mag_type'] = []
        d2['mag_sites'] = []
        d2['intFieldDB'] = []
        d2['intField'] = []
        d2['volume'] = []
        d2['magEl'] = []
        d2['magEl_a'] = []
        d2['magEl_b'] = []
        d2['magEl_c'] = []
                
        for (var i = 0; i < inds.length; i++) {
            d2['id'].push(d1['id'][inds[i]])
            d2['compound'].push(d1['compound'][inds[i]])
            d2['symmetry'].push(d1['symmetry'][inds[i]])
            d2['spacegroup'].push(d1['spacegroup'][inds[i]])
            d2['polymorphs'].push(d1['polymorphs'][inds[i]])
            d2['mag_type'].push(d1['mag_type'][inds[i]])
            d2['mag_sites'].push(d1['mag_sites'][inds[i]])
            d2['intFieldDB'].push(d1['intFieldDB'][inds[i]])
            d2['intField'].push(d1['intField'][inds[i]])
            d2['volume'].push(d1['volume'][inds[i]])
            d2['magEl'].push(d1['magEl'][inds[i]])
            d2['magEl_a'].push(d1['magEl_a'][inds[i]])
            d2['magEl_b'].push(d1['magEl_b'][inds[i]])
            d2['magEl_c'].push(d1['magEl_c'][inds[i]])            
        }
        s2.change.emit();
    """)
)

button.js_on_click(CustomJS(args=dict(source=s2), code=open(join(dirname(__file__), "download.js")).read()))


def select_entries():
    symType_val = symType.value
    compName_val = compName.value.strip()
    elemIn_val = elemIn.value.replace(" ", "")
    elemOut_val = elemOut.value.replace(" ", "")

    selected = datalist[
        (datalist.Mel_full >= magEl.value) &
        (datalist.magF_u >= intField.value)
        # (datalist.Year >= min_year.value) &
        # (datalist.Year <= max_year.value) &
    ]
    if (symType_val != "Any"):
        selected = selected[selected.lattice_system.str.contains(symType_val)==True]
    if (compName_val != ""):
        selected = selected[selected.pretty_formula.str.contains(compName_val)==True]
    if (elemIn_val != ""):
        contains_elements = set(str(elemIn_val).split(','))
        selected = selected[selected.species.str.replace("'", "").str.strip("[").str.strip("]").str.replace(" ", "").str.split(';').map(contains_elements.issubset)]
    if (elemOut_val != ""):
        not_elements = set(str(elemOut_val).split(','))
        selected = selected[selected.species.str.contains('|'.join(not_elements))==False]
    return selected


def update():
    df = select_entries()
    x_name = axis_map[x_axis.value]
    y_name = axis_map[y_axis.value]

    p.xaxis.axis_label = x_axis.value
    p.yaxis.axis_label = y_axis.value
    p.title.text = "%d compounds satisfy conditions" % len(df)
    source.data = dict(
        x=df[x_name],
        y=df[y_name],
        color=df["color"],
        id=df["ID"],
        compound=df["pretty_formula"],
        symmetry=df["lattice_system"],
        spacegroup=df["spacegroup"],
        polymorphs=df["polymorphs"],
        mag_type=df["mag_type"],
        mag_sites=df["mag_sites"],
        intFieldDB=df["mag_field"],
        intField=df["magF_u"],
        volume=df["volume_cell"],
        magEl=df["Mel_full"],
        magEl_a=df["Mel_a"],
        magEl_b=df["Mel_b"],
        magEl_c=df["Mel_c"],
    )
    source_back.data = dict(
        x=datalist[x_name],
        y=datalist[y_name],
        compound=datalist["pretty_formula"],
        symmetry=datalist["lattice_system"],
        intField=datalist["magF_u"],
        magEl=datalist["Mel_full"],
    )


controls = [symType, magEl, intField, compName, elemIn, elemOut, x_axis, y_axis]

for control in controls:
    control.on_change('value', lambda attr, old, new: update())

inputs = column(*controls, width=320, height=400, max_height=600)
inputs.sizing_mode = "fixed"
l = layout([[desc],
            [inputs, p, text],
            [button],
            [data_table]],
           sizing_mode="scale_both", max_width=1400, max_height=500)

update()  # initial load of the data
curdoc().add_root(l)
curdoc().title = "MCES"
