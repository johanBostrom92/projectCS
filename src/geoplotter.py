import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

#read a csv file and use semi-colon as delimeters for first row. 
df = pd.read_csv('C:\\Users\\Johan Boström\\Desktop\\projectcs_data.csv', delimiter=';')
df.head()

#where all the magic happens
fig = px.scatter_geo(df, size="popu" ,lat="lat", lon="long", hover_name="city", animation_frame="month", color="agent")


fig.update_layout(
        title_text = 'Corona spreading 2020',
        showlegend = True,
        geo = dict(
            scope = 'europe',
            landcolor = 'rgb(217, 217, 217)',
        )
    )

fig.show()


#Leagcy code
""" 
fig = go.Figure()


fig.add_trace(go.Scattergeo(
    #locationmode = 'europe',
    lon = df['long'],
    lat = df['lat'],
    animation_frame=df['month'],
        marker = dict(
        size = df['popu']/scale,
        line_color='rgb(40,40,40)',
        line_width=0.5,
        sizemode = 'area'
        ) ))
"""