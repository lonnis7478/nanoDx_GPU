import pandas
import pandas as pd
from plotnine import *
import numpy as np


df = pd.read_csv("../scripts/df_data.csv", index_col=0)

# ---------------------------- PLOT THE TSNE ----------------------------

def plot_tsne(df):

    # Read data from a file into a DataFrame
    colorMap = pd.read_csv("../static/colorMap_Capper_et_al.txt", skip_blank_lines=True, delimiter="\t")

    # Convert DataFrame to a Pandas DataFrame (equivalent to a tibble)
    colorMap = pd.DataFrame(colorMap)

    # Group by 'group' column
    grouped = colorMap.groupby('group')

    # Sort within each group by 'methylationClass' column
    colorMap = colorMap.groupby('group').apply(lambda x: x.sort_values('methylationClass'))

    # Add a row with 'color' set to 'white' at the beginning of each group

    colorMap = colorMap.apply(lambda x: x.append(pandas.Series({'color': 'white'}), ignore_index=True))
    colorMap['colorLabel'] = colorMap.apply(lambda row: f"**{row['group']}**" if pd.isna(row['methylationClass']) else row['methylationClass'], axis=1)
    colorMap['color'].fillna('grey', inplace=True)
    colorMap.loc[colorMap['colorLabel'] == 'unknown', 'color'] = 'red'
    hexCol = colorMap['color']
    hexCol.index = colorMap['colorLabel']

    colorMap["colorLabel"] = colorMap['colorLabel'].append(pandas.Series({"group":"unknown"}))

    df["Dx"].fillna("unknown", inplace=True)
    df['Dx'] = pd.Categorical(df['Dx'], categories=colorMap['colorLabel'])


    # Create the plot
    p = (ggplot(df.sort_values('Dx'),
                aes(x='X1', y='X2', color='Dx')) +
         geom_point() +
         coord_fixed(ratio=1) +
         theme_classic() +
         ggtitle("t-SNE, perplexity = 30") +
         scale_color_manual(values=hexCol) +
         scale_shape_manual(values=[16, 3]) +
         scale_size_manual(values=[1, 4]) +
         guides(colour=guide_legend(title="Methylation class",
                                    title_position="top",

                                    )) +
         theme(legend_text=element_text(size=4))
    )

    ggsave(p, filename='python_plot.png', dpi=500)
