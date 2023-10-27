from os.path import dirname, join
import numpy as np
import pandas as pd
from datetime import datetime
from datetime import date
from os import listdir
import sys
from bokeh.io import curdoc
from bokeh.document import Document
from bokeh.layouts import column, layout
from bokeh.models import ColumnDataSource,  DataTable, TableColumn, NumberFormatter, CustomJS, Panel, Tabs
from bokeh.models import Div, Select, Slider, TextInput, Button, Spinner, CheckboxGroup, CheckboxButtonGroup
from bokeh.models import HoverTool, CrosshairTool, HTMLTemplateFormatter
from bokeh.plotting import figure


def get_datasource(datafile):
    # datafile = 'outdir_out.csv'
    datalist = pd.read_csv(datafile)
    print("Reading data from ", datafile)
    datalist.reset_index(inplace=True)  # setting index to be a separate column
    datalist.fillna(0, inplace=True)  # replace missing values with zero
    datalist["color"] = np.where(datalist["Mel_full"] >= 0, "orange", "grey")  # adding a column with colour
    return datalist


dat_files = [f for f in listdir('.') if f.endswith('.csv')]
select_file = Select(title="Database file:", value=dat_files[0], options=dat_files)

# if len(dat_files) != 1:
#     raise ValueError('should be only one .csv file in the current directory')

datalist1 = get_datasource(dat_files[0])
# datalist2 = get_datasource(dat_files[1])
# datalist3 = get_datasource(dat_files[2])


def select_datasource():
    if select_file.value == dat_files[0]:
        datalist = datalist1
    # if select_file.value == dat_files[1]:
    #     datalist = datalist2
    # if select_file.value == dat_files[2]:
    #     datalist = datalist3
    back = (datalist.sort_values(by=axis_map[y_axis.value], ascending=True))  # making a separate dataframe for non-highlighted entries
    return datalist, back


# datalist = pd.read_csv(datafile)
# print("Reading data from ", datafile)
# datalist.reset_index(inplace=True)  # setting index to be a separate column
# datalist.fillna(0, inplace=True)   # replace missing values with zero
# datalist["color"] = np.where(datalist["Mel_full"] >= 0, "orange", "grey")  # adding a column with colour
# datalist['date_complete'] = pd.to_datetime(datalist['date_complete'], format='%d-%m-%y')
# datalist['date_complete'] = datalist['date_complete'].dt.date
# datalist["color"] = np.where(datalist["Number_of_def_points"] > 2, "red", "grey")  # adding a column with colour

# selecting available axis
axis_map = {
    "Magneto Elastic Full Uniaxial": "Mel_full",
    "Magneto Elastic Volumetric": "Mel_V",
    "Applied Field Response": "BEXT",
    "Magneto Elastic Uniaxial a": "Mel_a",
    "Magneto Elastic Uniaxial b": "Mel_b",
    "Magneto Elastic Uniaxial c": "Mel_c",
    "Internal Field Calculated": "magF_u",
    "Internal Field Maximal": "magF_max",
   # "Internal Field Database": "mag_field",
    "Volume": "volume_cell",
    # "Rare": "rare_content",
    "Price Estimate": "price"
    # "Index": "index"
}

# textboxes with content marked in html, separate files

# template = """<span href="#" data-toggle="tooltip" title="<%= value %>"><%= value %></span>"""HTMLTemplateFormatter(template=template)

desc = Div(text=open(join(dirname(__file__), "Title.html")).read())
text = Div(text=open(join(dirname(__file__), "Description.html")).read(), sizing_mode="stretch_height", margin=(0, 0, 0, 100), width=600, align='center')

symType_title = Div(text="""Lattice type(s):""", margin=[5, 0, 0, 45], width=250, sizing_mode="stretch_height")

# Create Input controls
# symType = Select(title="Symmetry", value="Any", options=open(join(dirname(__file__), 'symmetries.txt')).read().split())
symmetries = ['cubic', 'hexagonal', 'tetragonal', 'orthorhombic', 'monoclinic', 'rhombohedral', 'triclinic']
symType = CheckboxGroup(labels=symmetries, active=[], margin=[0, 0, 0, 5])
databases = ['ICSD', 'COD', 'MP']
databaseType = CheckboxButtonGroup(labels=databases, active=[0, 1, 2], margin=[0, 0, 0, 5])
magEl = Slider(title="Magneto Elastic Full Uniaxial", value=0, start=0, end=datalist1['Mel_full'].max(), step=0.5)
magEl_V = Slider(title="Magneto Elastic Volumetric", value=0, start=datalist1['Mel_V'].min(), end=datalist1['Mel_V'].max(), step=0.5)
# magEl = Slider(title="Magneto Elastic Full Uniaxial", value=0, start=0, end=100, step=0.5)
# magEl_V = Slider(title="Magneto Elastic Volumetric", value=0, start=0, end=100, step=0.5)
intField = Slider(title="Internal Field", start=0, end=3, value=0, step=0.05)
maxField = Slider(title="Internal Field (max)", start=0, end=3, value=0, step=0.05)
#appField = Slider(title="Applied Field Response", start=-2, end=2, value=0, step=0.01)
rareContent = Slider(title="Rare Elements allowed %", start=0, end=100, value=100, step=1)
magActive = Slider(title="Magnetically active atoms %", start=0, end=100, value=50, step=1)
prices = Spinner(title="Maximum allowed price (EUR/kg)", low=0, high=datalist1['price'].max(), step=0.5, value=datalist1['price'].max())
fitMel = Spinner(title="Mag.El. QOF", low=0, high=1, step=0.01, value=0)
fitE = Spinner(title="Energy QOF", low=0, high=1, step=0.01, value=0)

# prices = Spinner(title="Price Estimate Max EUR/kg", low=0, high=55000, step=0.5, value=55000)

# Slider(title="Price Estimate Max EUR/kg", start=0, end=datalist['price'].max(), value=datalist['price'].max(), step=1)
# dateComplete = DateRangeSlider(title="Date Calculated", start=datalist['date_complete'].min(), end=datalist['date_complete'].max(), value=(datalist['date_complete'].min(), datalist['date_complete'].max()), step=1)
id_is = TextInput(title="Database ID")

compName = TextInput(title="Compound name contains")
elemIn = TextInput(title="Contains Element(s)")
elemOut = TextInput(title="Does not contain Element(s)")
x_axis = Select(title="X Axis", options=sorted(axis_map.keys()), value="Magneto Elastic Full Uniaxial")
y_axis = Select(title="Y Axis", options=sorted(axis_map.keys()), value="Internal Field Calculated")
# load_button = Button(label="Load Data", button_type="default", height=50, width=150, sizing_mode="fixed", margin=(20, 40, 0, 100))
save_button = Button(label="Export Data", button_type="success", height=40, width=120, sizing_mode="fixed", align='center')
# stop_button = Button(label="Stop Server", button_type="danger", height=50, width=150, sizing_mode="fixed", margin=(20, 40, 0, 100))
# switching_slider = Slider(title="Cutoff", start=0, end=datalist[axis_map[x_axis.value]].max(), value=0, step=((datalist[axis_map[x_axis.value]].max())/100))


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
                                    mag_active=[],
                       #             mag_field=[],
                                    magF_u=[],
                                    magF_max=[],
                                    Mel_full=[],
                                    Mel_V=[],
                                    Mel_a=[],
                                    Mel_b=[],
                                    Mel_c=[],
                             #       BEXT=[],
                                    volume_cell=[],
                                    lattice_system=[],
                                    database=[],
                                    # date_complete=[],
                                    rare_content=[],
                                    price=[],
                                    Mel_worst_fit=[],
                                    E_worst_fit=[])
                          )
# Create second Data Source that will be used by the plot
source_back = ColumnDataSource(data=dict(x=[],
                                         y=[],
                                         # index=[],
                                         pretty_formula=[],
                                         magF_u=[],
                                         magF_max=[],
                                         Mel_full=[],
                                         Mel_V=[],
                                   #      BEXT=[],
                                         lattice_system=[],
                                         database=[])
                               )

# setting up additional tools
TOOLTIPS = [
    ("Compound", "@pretty_formula"),
    ("Symmetry", "@lattice_system"),
    ("Mag.el", "@Mel_full"),
    ("Mag.el V", "@Mel_V"),
    ("Internal field", "@magF_u"),
  #  ("A.F. Response", "@BEXT"),
    ("Database", "@database")
]

hover = HoverTool(names=['needshover'], tooltips=TOOLTIPS)
crosshair = CrosshairTool(toggleable=False)

# creating plot area
p = figure(plot_height=500, plot_width=500, title="", toolbar_location="right", tools="pan,wheel_zoom,box_select,box_zoom,lasso_select,tap,reset", active_scroll="wheel_zoom", tooltips=TOOLTIPS,
           sizing_mode="fixed", margin=(0, 0, 0, 50))

# adding extra tools to the plot
p.add_tools(crosshair)
p.add_tools(hover)
p.toolbar.active_inspect = [hover, crosshair]

# scatter plots
p.circle(x="x", y="y", source=source_back, size=4, color="black", line_color=None, fill_alpha=0.5)
p.circle(x="x", y="y", source=source, size=5, color="color", line_color=None, fill_alpha=0.9, name='needshover')
# formatting axis style
p.axis.axis_label_text_font_size = "16pt"
p.title.text_font_size = "16pt"
p.xaxis.major_label_text_font_size = "16pt"
p.yaxis.major_label_text_font_size = "16pt"

# p2 = figure(plot_height=450, plot_width=600, title="", toolbar_location="right", tools="pan,wheel_zoom,box_select,box_zoom,lasso_select,tap,reset,save", active_scroll="wheel_zoom", tooltips=TOOLTIPS, sizing_mode="fixed")
# p2.add_tools(crosshair)
# p2.add_tools(hover)
# p2.toolbar.active_inspect = [hover, crosshair]
# p2.circle(x="index", y="y", source=source_back, size=4, color="black", line_color=None, fill_alpha=0.5)
# # formatting axis style
# p2.axis.axis_label_text_font_size = "20pt"
# p2.title.text_font_size = "16pt"
# p2.xaxis.major_label_text_font_size = "20pt"
# p2.yaxis.major_label_text_font_size = "20pt"


# important_time = Span(location=2, dimension='width', line_color='red', line_dash='dashed', line_width=3)
# p2.add_layout(important_time)

# data source for table
s2 = ColumnDataSource(data=dict(ID=[],
                                pretty_formula=[],
                                spacegroup=[],
                                polymorphs=[],
                                mag_sites=[],
                                mag_active=[],
                            #    mag_field=[],
                                magF_u=[],
                                magF_max=[],
                                Mel_full=[],
                                Mel_V=[],
                                # Mel_a=[],
                                # Mel_b=[],
                                # Mel_c=[],
                              #  BEXT=[],
                                volume_cell=[],
                                lattice_system=[],
                                rare_content=[],
                                price=[],
                                database=[])
                      )
# table styling
columns = [
    TableColumn(field="ID", title="ID"),
    TableColumn(field="pretty_formula", title="Compound"),
    TableColumn(field="lattice_system", title="Symmetry"),
    TableColumn(field="spacegroup", title="Space group"),
    TableColumn(field="polymorphs", title="Polymorphs"),
    TableColumn(field="mag_sites", title="Sites"),
    TableColumn(field="mag_active", title="Mag.active", formatter=NumberFormatter(format="0.00")),
  #  TableColumn(field="mag_field", title="Internal Field DB", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="magF_u", title="Internal Field", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="magF_max", title="Internal Field (max)", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="volume_cell", title="Volume", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="Mel_full", title="Magneto Elastic", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="Mel_V", title="Mag.El. Volumetric", formatter=NumberFormatter(format="0.00")),
  #  TableColumn(field="BEXT", title="A.F. Response", formatter=NumberFormatter(format="0.00")),
    # TableColumn(field="rare_content", title="RE", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="price", title="Price", formatter=NumberFormatter(format="0.00")),
    TableColumn(field="database", title="Database"),

]
# creating table
data_table = DataTable(source=s2, columns=columns, sizing_mode="stretch_height", width=1100)

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
            d2['mag_sites'] = []
            d2['mag_active'] = []

            d2['magF_u'] = []
            d2['magF_max'] = []
            d2['volume_cell'] = []
            d2['Mel_full'] = []
            d2['Mel_V'] = []
            d2['Mel_a'] = []
            d2['Mel_b'] = []
            d2['Mel_c'] = []

            d2['rare_content'] = []
            d2['price'] = []
            d2['Mel_worst_fit'] = []
            d2['E_worst_fit'] = []
            d2['database'] = []

            for (var i = 0; i < inds.length; i++) {
                d2['ID'].push(d1['ID'][inds[i]])
                d2['pretty_formula'].push(d1['pretty_formula'][inds[i]])
                d2['lattice_system'].push(d1['lattice_system'][inds[i]])
                d2['spacegroup'].push(d1['spacegroup'][inds[i]])
                d2['polymorphs'].push(d1['polymorphs'][inds[i]])
                d2['mag_sites'].push(d1['mag_sites'][inds[i]])
                d2['mag_active'].push(d1['mag_active'][inds[i]])

                d2['magF_u'].push(d1['magF_u'][inds[i]])
                d2['magF_max'].push(d1['magF_max'][inds[i]])
                d2['volume_cell'].push(d1['volume_cell'][inds[i]])
                d2['Mel_full'].push(d1['Mel_full'][inds[i]])
                d2['Mel_V'].push(d1['Mel_V'][inds[i]])                
                d2['Mel_a'].push(d1['Mel_a'][inds[i]])
                d2['Mel_b'].push(d1['Mel_b'][inds[i]])
                d2['Mel_c'].push(d1['Mel_c'][inds[i]])

                d2['rare_content'].push(d1['rare_content'][inds[i]])            
                d2['price'].push(d1['price'][inds[i]])      
                d2['Mel_worst_fit'].push(d1['Mel_worst_fit'][inds[i]])    
                d2['E_worst_fit'].push(d1['E_worst_fit'][inds[i]])        
                d2['database'].push(d1['database'][inds[i]])

            }
            s2.change.emit();
        """
                                                 )
                             )


# populating the datasourse by values according to dynamic menu
def select_entries():
    # symType_val = symType.value
    datalist, back = select_datasource()
    symType_val = [symmetries[x] for x in symType.active]
    if len(symType_val) == 0:
        symType_val = "Any"
    databaseType_val = [databases[x] for x in databaseType.active]
    if len(databaseType_val) == 0:
        databaseType_val = "Any"
    compName_val = compName.value.strip()
    elemIn_val = elemIn.value.replace(" ", "")
    elemOut_val = elemOut.value.replace(" ", "")
    id_is_val = id_is.value.strip()

    selected = datalist[
        (datalist.Mel_full >= magEl.value) &
        (datalist.Mel_V >= magEl_V.value) &
        (datalist.magF_u >= intField.value) &
        (datalist.magF_max >= maxField.value) &
    #    (datalist.BEXT >= appField.value) &
        (datalist.rare_content <= rareContent.value) &
        (datalist.mag_active >= magActive.value) &
        (datalist.price <= prices.value) &
        (datalist.Mel_worst_fit >= fitMel.value) &
        (datalist.E_worst_fit >= fitE.value)  # &

        # (datalist['date_complete'] >= date.fromtimestamp(date_complete_val[0]/1000)) &
        # (datalist['date_complete'] <= date.fromtimestamp(date_complete_val[1]/1000))
        # (date.strptime(datalist.date_complete) >= dateComplete.value[0]) &
        # (datalist.date_complete <= dateComplete.value[1])
        # (datalist.Year >= min_year.value) &
        # (datalist.Year <= max_year.value) &
        ]
    if (symType_val != "Any"):
        # selected = selected[selected.lattice_system.str.contains(symType_val)==True]
        selected = selected[selected.lattice_system.isin(symType_val)]
    if (databaseType_val != "Any"):
        selected = selected[selected.database.isin(databaseType_val)]
    if (compName_val != ""):
        contains_compounds = set(str(compName_val).split(','))
        selected = selected[selected.pretty_formula.isin(contains_compounds)]
    if (id_is_val != ""):
        contains_ids = set(str(id_is_val).replace(" ", "").split(','))
        selected = selected[selected.ID.isin(contains_ids)]
    if (elemIn_val != ""):
        contains_elements = set(str(elemIn_val).split(','))
        selected = selected[selected.species.str.replace("'", "").str.strip("[").str.strip("]").str.replace(" ", "").str.split(';').map(contains_elements.issubset)]
    if (elemOut_val != ""):
        not_elements = set(str(elemOut_val).split(','))
        selected = selected[selected.species.str.contains('|'.join(not_elements)) == False]
    return selected, back


# update function for both plot datasources
def update():
    df, back = select_entries()
    x_name = axis_map[x_axis.value]
    y_name = axis_map[y_axis.value]

    p.xaxis.axis_label = x_axis.value
    p.yaxis.axis_label = y_axis.value
    # p.title.text = "%d compounds satisfy conditions \n(%d%%)" % (len(df), (100*(len(df)/len(datalist))))
    p.title.text = "%d compounds satisfy conditions" % len(df)

    # back2 = ((back.sort_values(by="index", ascending=True)))
    # back["index"] = back2.index.values

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
        mag_sites=df["mag_sites"],
        mag_active=df["mag_active"],
       #mag_field=df["mag_field"],
        magF_u=df["magF_u"],
        magF_max=df["magF_max"],
        volume_cell=df["volume_u"],
        Mel_full=df["Mel_full"],
        Mel_V=df["Mel_V"],
        Mel_a=df["Mel_a"],
        Mel_b=df["Mel_b"],
        Mel_c=df["Mel_c"],
   #     BEXT=df["BEXT"],
        # date_complete=df["date_complete"],
        rare_content=df["rare_content"],
        price=df["price"],
        Mel_worst_fit=df["Mel_worst_fit"],
        E_worst_fit=df["E_worst_fit"],
        database=df["database"]
    )
    source_back.data = dict(
        x=back[x_name],
        y=back[y_name],
        index=back["index"],
        pretty_formula=back["pretty_formula"],
        lattice_system=back["lattice_system"],
        magF_u=back["magF_u"],
        magF_max=back["magF_max"],
        Mel_full=back["Mel_full"],
        Mel_V=back["Mel_V"],
  #      BEXT=back["BEXT"],
        database=back["database"]
    )


# button event to save data, links to a custom js script
save_button.js_on_click(CustomJS(args=dict(source=source), code=open(join(dirname(__file__), "download.js")).read()))

# dynamic controls
controls1 = [select_file, x_axis, y_axis, id_is, compName, elemIn, elemOut, prices, databaseType]
for control1 in controls1:
    try:
        control1.on_change('value', lambda attr, old, new: update())
    except:
        databaseType.on_change('active', lambda attr, old, new: update())

controls2 = [symType, magEl, magEl_V, intField, maxField, magActive, rareContent]
for control2 in controls2:
    try:
        control2.on_change('value_throttled', lambda attr, old, new: update())
    except:
        symType.on_change('active', lambda attr, old, new: update())
#
controls3 = [fitMel, fitE]
for control3 in controls3:
    control3.on_change('value', lambda attr, old, new: update())

# creating input panels
inputs1 = column(*controls1, width=250)
inputs1.sizing_mode = "stretch_height"

inputs2 = column(*controls2, width=250, margin=(0, 0, 0, 40))
inputs2.sizing_mode = "stretch_height"

inputs3 = column(*controls3, width=250)
inputs3.sizing_mode = "fixed"

# created layouts
Layout1 = layout([
    [inputs1, p, [symType_title, inputs2, save_button], inputs3],
    [data_table]],
    sizing_mode="scale_both")
#
# Layout2 = layout([
#             [inputs3, [p2, save_button]]],
#                  sizing_mode="scale_both")

Layout3 = layout([text],
                 sizing_mode="scale_both")

tab1 = Panel(child=Layout1, title="Evaluation Plots")
# tab2 = Panel(child=Layout2, title="Sorted Plots")
tab3 = Panel(child=Layout3, title="Help")
tabs = Tabs(tabs=[tab1, tab3])

Layout0 = layout([
    [desc], [tabs]],
    sizing_mode="scale_both")

update()  # initial load of the data
curdoc().add_root(Layout0)
curdoc().title = "MCES"  # web page title