import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


#ref: https://plotly.com/python/bubble-maps/
#ref: https://plotly.com/python-api-reference/generated/plotly.express.scatter_geo.html


#read a csv file and use semi-colon as delimeters for first row.
df = pd.read_csv("..\\lib\\built_covid_data\\coviddata.csv", delimiter=';')
#df = pd.read_csv("C:\\Users\\Johan Boström\\Documents\GitHub\\projectCS\\lib\\built_covid_data\\coviddata.csv", delimiter=';')
df.head()


#create the actual map
fig = px.scatter_geo(df, size="popu" ,lat="lat", lon="long", hover_name="city", animation_frame="month", color="agent")


fig.update_layout(
        title_text = 'Corona spread 20XX',
        showlegend = True,
        geo = dict(
            scope = 'europe',
            landcolor = 'rgb(217, 217, 217)',
        )
    )

fig.show()
