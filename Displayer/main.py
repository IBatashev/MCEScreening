from os.path import dirname, join
import numpy as np
import pandas as pd

from bokeh.io import curdoc
from bokeh.layouts import column, layout
from bokeh.models import ColumnDataSource,  DataTable, TableColumn, NumberFormatter, CustomJS
from bokeh.models import Div, Select, Slider, TextInput, Button
from bokeh.models import HoverTool, CrosshairTool
from bokeh.plotting import figure

from tornado.ioloop import IOLoop

datalist = pd.read_csv('D:/MCES/MP/step1_success_sieved_out.csv')
datalist.reset_index(inplace=True) # setting index to be a separate column
datalist.fillna(0, inplace=True)   # replace missing values with zero

datalist["color"] = np.where(datalist["Mel_full"] >= 0, "orange", "grey") # adding a column with colour

# selecting available axis
axis_map = {
    "Magneto Elastic": "Mel_full",
    "Mag. El. a": "Mel_a",
    "Mag. El. b": "Mel_b",
    "Mag. El. c": "Mel_c",
    "Internal Field Calculated": "magF_u",
    "Internal Field Database": "mag_field",
    "Volume": "volume_cell",
    # "Index": "index"
}

# textboxes with content marked in html, separate files
desc = Div(text=open(join(dirname(__file__), "Title.html")).read(), sizing_mode="stretch_width", max_width=1400)
text = Div(text=open(join(dirname(__file__), "Description.html")).read(), sizing_mode="stretch_width", max_width=320, margin=(40, 40, 40, 40))

# Create Input controls
symType = Select(title="Symmetry", value="Any", options=open(join(dirname(__file__), 'symmetries.txt')).read().split())
magEl = Slider(title="Magneto Elastic", value=0, start=0, end=datalist['Mel_full'].max(), step=0.5)
intField = Slider(title="Internal Field", start=0, end=datalist['magF_u'].max(), value=0, step=0.05)
compName = TextInput(title="Compound name contains")
elemIn = TextInput(title="Contains Element(s)")
elemOut = TextInput(title="Does not contain Element(s)")
x_axis = Select(title="X Axis", options=sorted(axis_map.keys()), value="Magneto Elastic")
y_axis = Select(title="Y Axis", options=sorted(axis_map.keys()), value="Internal Field Calculated")
save_button = Button(label="Save Table", button_type="success", height=50, width=150, sizing_mode="fixed")
stop_button = Button(label="STOP", button_type="danger", height=40, width=75, sizing_mode="fixed", margin=(0, 40, 0, 175))

# Create Column Data Source that will be used by the plot
source = ColumnDataSource(data=dict(x=[],
                                    y=[],
                                    index=[],
                                    color=[],
                                    ID=[],
                                    pretty_formula=[],
                                    spacegroup=[],
                                    polymorphs=[],
                                    mag_sites=[],
                                    mag_type=[],
                                    mag_field=[],
                                    magF_u=[],
                                    Mel_full=[],
                                    Mel_a=[],
                                    Mel_b=[],
                                    Mel_c=[],
                                    volume_cell=[],
                                    lattice_system=[])
                          )
# Create second Data Source that will be used by the plot
source_back = ColumnDataSource(data=dict(x=[],
                                         y=[],
                                         pretty_formula=[],
                                         magF_u=[],
                                         Mel_full=[],
                                         lattice_system=[])
                               )

# setting up additional tools
TOOLTIPS = [
    ("Compound", "@pretty_formula"),
    ("Symmetry", "@lattice_system"),
    ("Mel", "@Mel_full"),
    ("Field", "@magF_u")
]

hover = HoverTool(names=['needshover'], tooltips=TOOLTIPS)
crosshair = CrosshairTool(toggleable=False)

# creating plot area
p = figure(plot_height=550, plot_width=550, title="", toolbar_location="right", tools="pan,wheel_zoom,box_select,box_zoom,lasso_select,tap,reset,save", active_scroll="wheel_zoom", tooltips=TOOLTIPS, sizing_mode="fixed")

# adding extra tools to the plot
p.add_tools(crosshair)
p.add_tools(hover)
p.toolbar.active_inspect = [hover, crosshair]

# scatter plots
p.circle(x="x", y="y", source=source_back, size=4, color="black", line_color=None, fill_alpha=0.5)
p.circle(x="x", y="y", source=source, size=5, color="color", line_color=None, fill_alpha=0.9, name='needshover')

# formatting axis style
p.axis.axis_label_text_font_size = "20pt"
p.title.text_font_size = "20pt"
p.xaxis.major_label_text_font_size = "20pt"
p.yaxis.major_label_text_font_size = "20pt"

# data source for table
s2 = ColumnDataSource(data=dict(ID=[],
                                pretty_formula=[],
                                spacegroup=[],
                                polymorphs=[],
                                mag_sites=[],
                                mag_type=[],
                                mag_field=[],
                                magF_u=[],
                                Mel_full=[],
                                Mel_a=[],
                                Mel_b=[],
                                Mel_c=[],
                                volume_cell=[],
                                lattice_system=[])
                      )
# table styling
columns = [
    TableColumn(field="ID", title="ID"),
    TableColumn(field="pretty_formula", title="Compound"),
    TableColumn(field="lattice_system", title="Symmetry"),
    TableColumn(field="spacegroup", title="Space group"),
    TableColumn(field="polymorphs", title="Polymorphs"),
    TableColumn(field="mag_type", title="Ordering"),
    TableColumn(field="mag_sites", title="Sites"),
    TableColumn(field="mag_field", title="Internal Field DB", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="magF_u", title="Internal Field", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="volume_cell", title="Volume", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="Mel_full", title="Magneto Elastic", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="Mel_a", title="Mag.El. a", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="Mel_b", title="Mag.El. b", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="Mel_c", title="Mag.El. c", formatter=NumberFormatter(format="0.00")),
]
# creating table
data_table = DataTable(source=s2, columns=columns, width=600)

# actions for selection
source.selected.js_on_change('indices', CustomJS(args=dict(s1=source, s2=s2), code="""
        var inds = cb_obj.indices;
        var d1 = s1.data;
        var d2 = s2.data;
        d2['ID'] = []
        d2['pretty_formula'] = []
        d2['lattice_system'] = []
        d2['spacegroup'] = []
        d2['polymorphs'] = []
        d2['mag_type'] = []
        d2['mag_sites'] = []
        d2['mag_field'] = []
        d2['magF_u'] = []
        d2['volume_cell'] = []
        d2['Mel_full'] = []
        d2['Mel_a'] = []
        d2['Mel_b'] = []
        d2['Mel_c'] = []
                
        for (var i = 0; i < inds.length; i++) {
            d2['ID'].push(d1['ID'][inds[i]])
            d2['pretty_formula'].push(d1['pretty_formula'][inds[i]])
            d2['lattice_system'].push(d1['lattice_system'][inds[i]])
            d2['spacegroup'].push(d1['spacegroup'][inds[i]])
            d2['polymorphs'].push(d1['polymorphs'][inds[i]])
            d2['mag_type'].push(d1['mag_type'][inds[i]])
            d2['mag_sites'].push(d1['mag_sites'][inds[i]])
            d2['mag_field'].push(d1['mag_field'][inds[i]])
            d2['magF_u'].push(d1['magF_u'][inds[i]])
            d2['volume_cell'].push(d1['volume_cell'][inds[i]])
            d2['Mel_full'].push(d1['Mel_full'][inds[i]])
            d2['Mel_a'].push(d1['Mel_a'][inds[i]])
            d2['Mel_b'].push(d1['Mel_b'][inds[i]])
            d2['Mel_c'].push(d1['Mel_c'][inds[i]])            
        }
        s2.change.emit();
    """
                                                 )
                             )


# populating the datasourse by values according to dynamic menu
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


# update function for both plot datasources
def update():
    df = select_entries()
    x_name = axis_map[x_axis.value]
    y_name = axis_map[y_axis.value]

    p.xaxis.axis_label = x_axis.value
    p.yaxis.axis_label = y_axis.value
    p.title.text = "%d compounds satisfy conditions" % len(df)

    # df2 = ((df.sort_values(by=y_name, ascending=True).reset_index()))
    # df3 = ((df2.sort_values(by="index", ascending=True)))
    # df["index"] = df3.index.values

    source.data = dict(
        x=df[x_name],
        y=df[y_name],
        # index=df["index"],
        color=df["color"],
        ID=df["ID"],
        pretty_formula=df["pretty_formula"],
        lattice_system=df["lattice_system"],
        spacegroup=df["spacegroup"],
        polymorphs=df["polymorphs"],
        mag_type=df["mag_type"],
        mag_sites=df["mag_sites"],
        mag_field=df["mag_field"],
        magF_u=df["magF_u"],
        volume_cell=df["volume_cell"],
        Mel_full=df["Mel_full"],
        Mel_a=df["Mel_a"],
        Mel_b=df["Mel_b"],
        Mel_c=df["Mel_c"],
    )
    source_back.data = dict(
        x=datalist[x_name],
        y=datalist[y_name],
        index=datalist["index"],

        pretty_formula=datalist["pretty_formula"],
        lattice_system=datalist["lattice_system"],
        magF_u=datalist["magF_u"],
        Mel_full=datalist["Mel_full"],
    )


# function to stop server to use in a button event
def stop_server():
    IOLoop.current().stop()


# button event to save data, links to a custom js script
save_button.js_on_click(CustomJS(args=dict(source=s2), code=open(join(dirname(__file__), "download.js")).read()))
# button for stopping the server
stop_button.on_click(stop_server)

# dynamic controls
controls = [symType, magEl, intField, compName, elemIn, elemOut, x_axis, y_axis]
for control in controls:
    control.on_change('value', lambda attr, old, new: update())

# creating input panel
inputs = column(*controls, width=320, height=400, max_height=600)
inputs.sizing_mode = "fixed"

# created layout
l = layout([[desc],
            [inputs, p, [text, stop_button]],
            [save_button],
            [data_table]],
           sizing_mode="scale_both", max_width=1400, max_height=500)


update()  # initial load of the data
curdoc().add_root(l)
curdoc().title = "MCES"  # web page title
