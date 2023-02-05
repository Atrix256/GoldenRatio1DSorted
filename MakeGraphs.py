import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

fileNames = ["Sine", "Step", "Triangle"]

for fileName in fileNames:
    df = pd.read_csv("out/" + fileName + ".csv").drop(['Samples'], axis=1)

    fig, ax = plt.subplots()
    for col in df.columns:
        ax.plot(df[col], label=col)

    ax.legend()

    plt.title(fileName)
    plt.ylabel('RMSE')
    plt.xlabel('Samples')

    fig.axes[0].set_xscale('log', base=2)
    fig.axes[0].set_yscale('log', base=2)

    fig.axes[0].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    fig.axes[0].get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    fig.tight_layout()
    fig.savefig("out/" + fileName + ".png", bbox_inches='tight')
