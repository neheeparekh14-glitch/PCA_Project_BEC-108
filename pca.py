import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import gzip


def get_gene_id(df_expr_cols, gene_name, filepath='columns.tsv.gz'):
    """Helper function to safely extract a gene ID from the mapping file."""
    valid_cols = set(df_expr_cols)
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        for line in f:
            if gene_name in line:
                parts = [p.strip(' "\n\r') for p in line.split('\t')]
                for part in parts:
                    if part in valid_cols and part.isdigit():
                        return part
    return None


def main():
    print("Loading data...")
    # 1. Load the labels and expression matrix
    df_class = pd.read_csv('class.tsv', sep='\t', header=None)
    labels = df_class.iloc[:, -1].values

    df_expr = pd.read_csv('filtered.tsv.gz', sep='\t')
    df_expr.columns = [str(col).strip() for col in df_expr.columns]  # Clean column spaces

    # 2. Extract Gene IDs
    xbp1_id = '4404'
    gata3_id = get_gene_id(df_expr.columns, 'GATA3')
    ccnb2_id = get_gene_id(df_expr.columns, 'CCNB2')  # Needed for Plot E

    if not gata3_id or not ccnb2_id:
        print("Warning: Could not find dynamic IDs, falling back to typical GEO IDs if applicable.")

    print(f"IDs - XBP1: {xbp1_id}, GATA3: {gata3_id}, CCNB2: {ccnb2_id}")

    # 3. Prepare Data Matrices
    # 2D Data (for plots a, b, c)
    X_2d = df_expr[[gata3_id, xbp1_id]].values
    gata3_expr = X_2d[:, 0]
    xbp1_expr = X_2d[:, 1]

    # Full Data (for plots d, e, f)
    X_full = df_expr.values

    # Masks for ER labels
    er_plus = (labels == 1)
    er_minus = (labels == 0)

    # 4. Run PCAs
    # PCA on 2 Genes
    pca_2d = PCA(n_components=2)
    X_2d_pca = pca_2d.fit_transform(X_2d)

    # PCA on All Genes
    pca_full = PCA()
    X_full_pca = pca_full.fit_transform(X_full)

    # ==========================================
    # PLOTTING FUNCTIONS
    # ==========================================

    # Setup the combined figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()

    def clean_axes(ax):
        """Removes top and right spines for the Nature aesthetic."""
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # --- Plot A: 2D Scatter ---
    print("Generating Plot A...")
    axes[0].scatter(gata3_expr[er_minus], xbp1_expr[er_minus], c='black', s=15, marker='s')
    axes[0].scatter(gata3_expr[er_plus], xbp1_expr[er_plus], c='red', s=15, marker='s')
    axes[0].set_xlabel('GATA3')
    axes[0].set_ylabel('XBP1')
    axes[0].set_title('a', loc='left', fontweight='bold', fontsize=16)
    clean_axes(axes[0])

    # Save individual
    fig_a = plt.figure(figsize=(5, 4));
    ax_a = fig_a.add_subplot(111)
    ax_a.scatter(gata3_expr[er_minus], xbp1_expr[er_minus], c='black', s=15, marker='s')
    ax_a.scatter(gata3_expr[er_plus], xbp1_expr[er_plus], c='red', s=15, marker='s')
    ax_a.set_xlabel('GATA3');
    ax_a.set_ylabel('XBP1')
    clean_axes(ax_a);
    fig_a.savefig('plot_a.png', dpi=300, bbox_inches='tight');
    plt.close(fig_a)

    # --- Plot B: 2D Scatter with PC Vectors ---
    print("Generating Plot B...")
    axes[1].scatter(gata3_expr[er_minus], xbp1_expr[er_minus], c='black', s=15, marker='s', zorder=1)
    axes[1].scatter(gata3_expr[er_plus], xbp1_expr[er_plus], c='red', s=15, marker='s', zorder=1)

    mean_x, mean_y = pca_2d.mean_

    # Dynamic scaling: Make the vector length relative to the data's actual visual spread
    spread_x = np.max(gata3_expr) - np.min(gata3_expr)
    spread_y = np.max(xbp1_expr) - np.min(xbp1_expr)
    max_spread = max(spread_x, spread_y)

    # Plotting PC lines reliably using plot() instead of finicky arrows
    for i, (comp, var_ratio) in enumerate(zip(pca_2d.components_, pca_2d.explained_variance_ratio_)):
        # Scale length by the spread of the data and how much variance it explains.
        # (Adding 0.3 ensures PC2 doesn't become microscopically short if PC1 dominates)
        vector_length = max_spread * (var_ratio + 0.3) * 0.7
        vector = comp * vector_length

        # Draw the solid line across the mean
        axes[1].plot([mean_x - vector[0], mean_x + vector[0]],
                     [mean_y - vector[1], mean_y + vector[1]],
                     color='black', lw=1.5, zorder=2)

        # Add the text label slightly past the end of the line
        axes[1].text(mean_x + vector[0] * 1.15, mean_y + vector[1] * 1.15, f'PC{i + 1}',
                     rotation=(np.degrees(np.arctan2(vector[1], vector[0]))),
                     ha='center', va='center', fontweight='bold', zorder=3)

    axes[1].set_xlabel('GATA3')
    axes[1].set_ylabel('XBP1')
    axes[1].set_title('b', loc='left', fontweight='bold', fontsize=16)
    clean_axes(axes[1])

    # (Don't forget to update the single saved figure code for plot_b.png as well)
    fig_b = plt.figure(figsize=(5, 4));
    ax_b = fig_b.add_subplot(111)
    ax_b.scatter(gata3_expr[er_minus], xbp1_expr[er_minus], c='black', s=15, marker='s', zorder=1)
    ax_b.scatter(gata3_expr[er_plus], xbp1_expr[er_plus], c='red', s=15, marker='s', zorder=1)
    for i, (comp, var_ratio) in enumerate(zip(pca_2d.components_, pca_2d.explained_variance_ratio_)):
        vector_length = max_spread * (var_ratio + 0.3) * 0.7
        vector = comp * vector_length
        ax_b.plot([mean_x - vector[0], mean_x + vector[0]], [mean_y - vector[1], mean_y + vector[1]], color='black',
                  lw=1.5, zorder=2)
        ax_b.text(mean_x + vector[0] * 1.15, mean_y + vector[1] * 1.15, f'PC{i + 1}',
                  rotation=(np.degrees(np.arctan2(vector[1], vector[0]))), ha='center', va='center', fontweight='bold')
    ax_b.set_xlabel('GATA3');
    ax_b.set_ylabel('XBP1');
    clean_axes(ax_b)
    fig_b.savefig('plot_b.png', dpi=300, bbox_inches='tight');
    plt.close(fig_b)

    # --- Plot C: 1D Projection ---
    print("Generating Plot C...")
    pc1_1d = X_2d_pca[:, 0]
    axes[2].axhline(3, color='gray', lw=1, zorder=1);
    axes[2].axhline(2, color='gray', lw=1, zorder=1);
    axes[2].axhline(1, color='gray', lw=1, zorder=1)

    axes[2].scatter(pc1_1d[er_minus], [3] * sum(er_minus), c='black', s=15, zorder=2)
    axes[2].scatter(pc1_1d[er_plus], [3] * sum(er_plus), c='red', s=15, zorder=2)
    axes[2].scatter(pc1_1d[er_minus], [2] * sum(er_minus), c='black', s=15, zorder=2)
    axes[2].scatter(pc1_1d[er_plus], [1] * sum(er_plus), c='red', s=15, zorder=2)

    axes[2].set_yticks([1, 2, 3]);
    axes[2].set_yticklabels(['ER+', 'ER-', 'All'])
    axes[2].set_xlabel('Projection onto PC1')
    axes[2].set_title('c', loc='left', fontweight='bold', fontsize=16)
    clean_axes(axes[2]);
    axes[2].spines['left'].set_visible(False);
    axes[2].tick_params(axis='y', length=0)

    fig_c = plt.figure(figsize=(5, 4));
    ax_c = fig_c.add_subplot(111)
    ax_c.axhline(3, color='gray', lw=1, zorder=1);
    ax_c.axhline(2, color='gray', lw=1, zorder=1);
    ax_c.axhline(1, color='gray', lw=1, zorder=1)
    ax_c.scatter(pc1_1d[er_minus], [3] * sum(er_minus), c='black', s=15, zorder=2);
    ax_c.scatter(pc1_1d[er_plus], [3] * sum(er_plus), c='red', s=15, zorder=2)
    ax_c.scatter(pc1_1d[er_minus], [2] * sum(er_minus), c='black', s=15, zorder=2);
    ax_c.scatter(pc1_1d[er_plus], [1] * sum(er_plus), c='red', s=15, zorder=2)
    ax_c.set_yticks([1, 2, 3]);
    ax_c.set_yticklabels(['ER+', 'ER-', 'All']);
    ax_c.set_xlabel('Projection onto PC1')
    clean_axes(ax_c);
    ax_c.spines['left'].set_visible(False);
    ax_c.tick_params(axis='y', length=0)
    fig_c.savefig('plot_c.png', dpi=300, bbox_inches='tight');
    plt.close(fig_c)

    # --- Plot D: Variance Explained (Scree Plot) ---
    print("Generating Plot D...")
    var_exp = pca_full.explained_variance_ratio_ * 100
    axes[3].bar(range(len(var_exp)), var_exp, color='gray')
    axes[3].set_xlabel('Principal component')
    axes[3].set_ylabel('Proportion of variance (%)')
    axes[3].set_title('d', loc='left', fontweight='bold', fontsize=16)
    clean_axes(axes[3])

    fig_d = plt.figure(figsize=(5, 4));
    ax_d = fig_d.add_subplot(111)
    ax_d.bar(range(len(var_exp)), var_exp, color='gray')
    ax_d.set_xlabel('Principal component');
    ax_d.set_ylabel('Proportion of variance (%)');
    clean_axes(ax_d)
    fig_d.savefig('plot_d.png', dpi=300, bbox_inches='tight');
    plt.close(fig_d)

    # --- Plot E: Full PCA Projection with Feature Loadings ---
    print("Generating Plot E...")
    pc1_full = X_full_pca[:, 0]
    pc2_full = X_full_pca[:, 1]

    axes[4].scatter(pc1_full[er_minus], pc2_full[er_minus], c='black', s=15, marker='s')
    axes[4].scatter(pc1_full[er_plus], pc2_full[er_plus], c='red', s=15, marker='s')
    axes[4].set_xlabel('Projection onto PC1')
    axes[4].set_ylabel('Projection onto PC2')
    axes[4].set_title('e', loc='left', fontweight='bold', fontsize=16)
    clean_axes(axes[4])

    # Adding gene loading vectors for XBP1 and CCNB2
    if ccnb2_id and gata3_id:
        idx_xbp1 = list(df_expr.columns).index(xbp1_id)
        idx_ccnb2 = list(df_expr.columns).index(ccnb2_id)

        # Extract loadings and scale for visibility
        loadings = pca_full.components_.T * np.sqrt(pca_full.explained_variance_) * 400

        for idx, name in zip([idx_xbp1, idx_ccnb2], ['XBP1', 'CCNB2']):
            axes[4].plot([0, loadings[idx, 0]], [0, loadings[idx, 1]], color='black', lw=1)
            axes[4].scatter(loadings[idx, 0], loadings[idx, 1], color='limegreen', zorder=3)
            axes[4].text(loadings[idx, 0] * 1.1, loadings[idx, 1] * 1.1, name)

    fig_e = plt.figure(figsize=(5, 4));
    ax_e = fig_e.add_subplot(111)
    ax_e.scatter(pc1_full[er_minus], pc2_full[er_minus], c='black', s=15, marker='s')
    ax_e.scatter(pc1_full[er_plus], pc2_full[er_plus], c='red', s=15, marker='s')
    if ccnb2_id and gata3_id:
        for idx, name in zip([idx_xbp1, idx_ccnb2], ['XBP1', 'CCNB2']):
            ax_e.plot([0, loadings[idx, 0]], [0, loadings[idx, 1]], color='black', lw=1)
            ax_e.scatter(loadings[idx, 0], loadings[idx, 1], color='limegreen', zorder=3)
            ax_e.text(loadings[idx, 0] * 1.1, loadings[idx, 1] * 1.1, name)
    ax_e.set_xlabel('Projection onto PC1');
    ax_e.set_ylabel('Projection onto PC2');
    clean_axes(ax_e)
    fig_e.savefig('plot_e.png', dpi=300, bbox_inches='tight');
    plt.close(fig_e)

    # --- Plot F: Full PCA with Proxy Subtypes (KMeans) ---
    print("Generating Plot F...")
    # Approximating missing Subtype labels using K-means
    kmeans = KMeans(n_clusters=3, random_state=42).fit(X_full)
    colors = {0: 'goldenrod', 1: 'royalblue', 2: 'limegreen'}
    c_labels = [colors[label] for label in kmeans.labels_]

    axes[5].scatter(pc1_full, pc2_full, c=c_labels, s=15, marker='s')
    axes[5].set_xlabel('Projection onto PC1')
    axes[5].set_ylabel('Projection onto PC2')
    axes[5].set_title('f', loc='left', fontweight='bold', fontsize=16)
    clean_axes(axes[5])

    fig_f = plt.figure(figsize=(5, 4));
    ax_f = fig_f.add_subplot(111)
    ax_f.scatter(pc1_full, pc2_full, c=c_labels, s=15, marker='s')
    ax_f.set_xlabel('Projection onto PC1');
    ax_f.set_ylabel('Projection onto PC2');
    clean_axes(ax_f)
    fig_f.savefig('plot_f.png', dpi=300, bbox_inches='tight');
    plt.close(fig_f)

    # --- Save and Show Combined Plot ---
    print("Saving combined plot...")
    plt.tight_layout()
    fig.savefig('plot_combined.png', dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    main()