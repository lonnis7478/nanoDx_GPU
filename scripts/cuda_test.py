import pandas as pd
from plotnine import *


df = pd.read_csv("data.csv", index_col=0)

# ---------------------------- PLOT THE TSNE ----------------------------

# Read data from a file into a DataFrame
colorMap = pd.read_csv("../static/colorMap_Capper_et_al.txt", skip_blank_lines=True, delimiter="\t")

# Convert DataFrame to a Pandas DataFrame (equivalent to a tibble)
colorMap = pd.DataFrame(colorMap)

# Group by 'group' column
grouped = colorMap.groupby('group')

# Sort within each group by 'methylationClass' column
colorMap = colorMap.groupby('group').apply(lambda x: x.sort_values('methylationClass'))

# Add a row with 'color' set to 'white' at the beginning of each group
colorMap = grouped.apply(lambda x: x.append({'color': 'white'}, ignore_index=True))

# Create the 'colorLabel' column based on 'methylationClass' and 'group'
colorMap['colorLabel'] = colorMap.apply(lambda row: f"**{row['group']}**" if pd.isna(row['methylationClass']) else row['methylationClass'], axis=1)

# Replace missing values in 'color' column with 'grey'
colorMap['color'].fillna('grey', inplace=True)

# Set 'unknown' color to 'red'
colorMap.loc[colorMap['colorLabel'] == 'unknown', 'color'] = 'red'

# Extract the 'color' column to 'hexCol'
hexCol = colorMap['color']

# Rename the elements in 'hexCol' with 'colorLabel' values
hexCol.index = colorMap['colorLabel']

df['Dx'] = pd.Categorical(df['Dx'], categories=colorMap['colorLabel'])


# Create the plot
p = (ggplot(df.query('Dx == "unknown"').sort_values('Dx'),
            aes(x='X1', y='X2', color='Dx', shape='Dx=="unknown", size=Dx=="unknown"')) +
     geom_point() +
     theme_classic() +
     ggtitle("t-SNE, perplexity = 30") +
     scale_color_manual(values=hexCol) +
     scale_shape_manual(values=[16, 3]) +
     scale_size_manual(values=[1, 4]) +
     guides(colour=guide_legend(title="Methylation class",
                                title_position="top",
                                override_aes={'shape': 15, 'size': 3}
                                )) +
     theme(legend_text=element_text(size=7))
)

ggsave(p, filename='python_plot.png', dpi=300)