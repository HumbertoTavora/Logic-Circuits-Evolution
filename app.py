import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- RAW DATA ---
# Each row corresponds to one genome (4PG1 ... 4PG9). For every input vector of the
# truth table we store a tuple: (input_index, input_bits, correct_responses), where
# correct_responses is the number of fault-injection trials that preserved the
# expected output for that input.
circuit_data = [
    [(1, [0, 0, 0, 0], 16), (2, [0, 0, 0, 1], 17), (3, [0, 0, 1, 0], 16), (4, [0, 0, 1, 1], 18), (5, [0, 1, 0, 0], 13), (6, [0, 1, 0, 1], 19), (7, [0, 1, 1, 0], 15), (8, [0, 1, 1, 1], 22), (9, [1, 0, 0, 0], 15), (10, [1, 0, 0, 1], 18), (11, [1, 0, 1, 0], 16), (12, [1, 0, 1, 1], 16), (13, [1, 1, 0, 0], 14), (14, [1, 1, 0, 1], 20), (15, [1, 1, 1, 0], 14), (16, [1, 1, 1, 1], 19)],
    [(1, [0, 0, 0, 0], 15), (2, [0, 0, 0, 1], 16), (3, [0, 0, 1, 0], 17), (4, [0, 0, 1, 1], 17), (5, [0, 1, 0, 0], 17), (6, [0, 1, 0, 1], 17), (7, [0, 1, 1, 0], 17), (8, [0, 1, 1, 1], 17), (9, [1, 0, 0, 0], 18), (10, [1, 0, 0, 1], 14), (11, [1, 0, 1, 0], 17), (12, [1, 0, 1, 1], 16), (13, [1, 1, 0, 0], 17), (14, [1, 1, 0, 1], 16), (15, [1, 1, 1, 0], 20), (16, [1, 1, 1, 1], 16)],
    [(1, [0, 0, 0, 0], 14), (2, [0, 0, 0, 1], 15), (3, [0, 0, 1, 0], 15), (4, [0, 0, 1, 1], 15), (5, [0, 1, 0, 0], 13), (6, [0, 1, 0, 1], 15), (7, [0, 1, 1, 0], 16), (8, [0, 1, 1, 1], 13), (9, [1, 0, 0, 0], 13), (10, [1, 0, 0, 1], 15), (11, [1, 0, 1, 0], 16), (12, [1, 0, 1, 1], 13), (13, [1, 1, 0, 0], 16), (14, [1, 1, 0, 1], 15), (15, [1, 1, 1, 0], 15), (16, [1, 1, 1, 1], 17)],
    [(1, [0, 0, 0, 0], 9), (2, [0, 0, 0, 1], 9), (3, [0, 0, 1, 0], 9), (4, [0, 0, 1, 1], 9), (5, [0, 1, 0, 0], 8), (6, [0, 1, 0, 1], 8), (7, [0, 1, 1, 0], 8), (8, [0, 1, 1, 1], 8), (9, [1, 0, 0, 0], 9), (10, [1, 0, 0, 1], 9), (11, [1, 0, 1, 0], 9), (12, [1, 0, 1, 1], 9), (13, [1, 1, 0, 0], 8), (14, [1, 1, 0, 1], 8), (15, [1, 1, 1, 0], 8), (16, [1, 1, 1, 1], 8)],
    [(1, [0, 0, 0, 0], 13), (2, [0, 0, 0, 1], 14), (3, [0, 0, 1, 0], 15), (4, [0, 0, 1, 1], 15), (5, [0, 1, 0, 0], 14), (6, [0, 1, 0, 1], 15), (7, [0, 1, 1, 0], 15), (8, [0, 1, 1, 1], 17), (9, [1, 0, 0, 0], 12), (10, [1, 0, 0, 1], 14), (11, [1, 0, 1, 0], 13), (12, [1, 0, 1, 1], 17), (13, [1, 1, 0, 0], 13), (14, [1, 1, 0, 1], 14), (15, [1, 1, 1, 0], 16), (16, [1, 1, 1, 1], 15)],
    [(1, [0, 0, 0, 0], 8), (2, [0, 0, 0, 1], 8), (3, [0, 0, 1, 0], 8), (4, [0, 0, 1, 1], 8), (5, [0, 1, 0, 0], 8), (6, [0, 1, 0, 1], 10), (7, [0, 1, 1, 0], 8), (8, [0, 1, 1, 1], 10), (9, [1, 0, 0, 0], 8), (10, [1, 0, 0, 1], 9), (11, [1, 0, 1, 0], 8), (12, [1, 0, 1, 1], 9), (13, [1, 1, 0, 0], 8), (14, [1, 1, 0, 1], 9), (15, [1, 1, 1, 0], 8), (16, [1, 1, 1, 1], 9)],
    [(1, [0, 0, 0, 0], 9), (2, [0, 0, 0, 1], 10), (3, [0, 0, 1, 0], 10), (4, [0, 0, 1, 1], 9), (5, [0, 1, 0, 0], 9), (6, [0, 1, 0, 1], 9), (7, [0, 1, 1, 0], 9), (8, [0, 1, 1, 1], 9), (9, [1, 0, 0, 0], 10), (10, [1, 0, 0, 1], 9), (11, [1, 0, 1, 0], 9), (12, [1, 0, 1, 1], 10), (13, [1, 1, 0, 0], 9), (14, [1, 1, 0, 1], 9), (15, [1, 1, 1, 0], 9), (16, [1, 1, 1, 1], 9)],
    [(1, [0, 0, 0, 0], 9), (2, [0, 0, 0, 1], 11), (3, [0, 0, 1, 0], 11), (4, [0, 0, 1, 1], 10), (5, [0, 1, 0, 0], 11), (6, [0, 1, 0, 1], 10), (7, [0, 1, 1, 0], 9), (8, [0, 1, 1, 1], 11), (9, [1, 0, 0, 0], 11), (10, [1, 0, 0, 1], 9), (11, [1, 0, 1, 0], 11), (12, [1, 0, 1, 1], 11), (13, [1, 1, 0, 0], 11), (14, [1, 1, 0, 1], 11), (15, [1, 1, 1, 0], 11), (16, [1, 1, 1, 1], 9)],
    [(1, [0, 0, 0, 0], 7), (2, [0, 0, 0, 1], 7), (3, [0, 0, 1, 0], 7), (4, [0, 0, 1, 1], 7), (5, [0, 1, 0, 0], 7), (6, [0, 1, 0, 1], 7), (7, [0, 1, 1, 0], 7), (8, [0, 1, 1, 1], 7), (9, [1, 0, 0, 0], 7), (10, [1, 0, 0, 1], 7), (11, [1, 0, 1, 0], 7), (12, [1, 0, 1, 1], 7), (13, [1, 1, 0, 0], 7), (14, [1, 1, 0, 1], 7), (15, [1, 1, 1, 0], 7), (16, [1, 1, 1, 1], 7)]
]
# Total number of logic gates in each genome (circuit size).
circuit_sizes = [23, 23, 21, 14, 20, 15, 15, 16, 12]

# --- PROCESSING ---
genome_labels = [f"4PG{i}" for i in range(1, 10)]
input_labels = [f"{i:04b}" for i in range(16)]

# Build matrices: rows = truth-table inputs, columns = genomes.
hits_matrix = np.array([[t[2] for t in c] for c in circuit_data]).T
annotation_matrix = np.array([[f"{h}/{circuit_sizes[idx]}" for h in [t[2] for t in c]]
                              for idx, c in enumerate(circuit_data)]).T

# Per-genome (column-wise) min-max normalisation to enable relative comparison
# of robustness across inputs within each circuit.
z_norm = np.zeros(hits_matrix.shape)
for col in range(hits_matrix.shape[1]):
    v_min, v_max = hits_matrix[:, col].min(), hits_matrix[:, col].max()
    if v_max > v_min:
        z_norm[:, col] = (hits_matrix[:, col] - v_min) / (v_max - v_min)
    else:
        z_norm[:, col] = 0.5

# --- FIGURE GENERATION (MATPLOTLIB) ---
def generate_eps():
    fig, ax = plt.subplots(figsize=(12, 10))

    # Draw the heatmap.
    sns.heatmap(
        z_norm,
        annot=annotation_matrix,
        fmt="",
        cmap="RdYlGn",
        xticklabels=genome_labels,
        yticklabels=input_labels,
        ax=ax,
        cbar_kws={'label': 'Relative robustness (column-normalised)'}
    )

    ax.set_title("Relative Robustness: Per-Genome Normalised Fault Tolerance",
                 fontsize=14, pad=20)
    ax.set_xlabel("Genome")
    ax.set_ylabel("Truth-table input (b₃b₂b₁b₀)")

    # Save as EPS. bbox_inches='tight' prevents labels from being clipped.
    plt.savefig("grafico_robustezz.eps", format='eps', bbox_inches='tight')
    return fig

# Streamlit interface
st.title("🛡️ EPS Generator: Circuit Robustness")

if st.button("🚀 Generate and Save EPS"):
    matplot_fig = generate_eps()
    st.pyplot(matplot_fig)
    st.success("✅ File 'grafico_robustez.eps' has been generated in the project folder!")

    # Provide the generated figure as a downloadable file.
    with open("grafico_robustez.eps", "rb") as f:
        st.download_button(
            label="📥 Download .EPS file",
            data=f,
            file_name="grafico_robustez.eps",
            mime="application/postscript"
        )
