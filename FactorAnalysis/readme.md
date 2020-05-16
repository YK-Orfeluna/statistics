# FactorAnalysis

## Enviornment

- Python (3.6.9)
- numpy (1.17.4)
- pandas (0.25.3)
- matplotlib (3.1.1)
- seaborn (0.9.0)
- factor_analyzer (0.3.2)

## Arguments

### Requires

- `inf`: Input file (csv or tsv only).
- `outd`: Output files will be saved in $outd.
- `N_adj`: Number of adjective pairs.

### Options

- `--inf_header`: If $inf has header, you have to activate this argument.
- `--inf_index`: If $inf has index, you have to activate this argument.
- `--inf_drop_NA`: If activate, NAs in $inf is removed.
- `--p_alpha`: Significance level for chi-squared test; default is 0.05.
- `--FA_rotation`: To choose rotation of factor analysis; default is None.
    - You can choose 'varimax', 'promax', 'oblimin', 'oblimax', 'quartimin', 'quartimax' or 'equamax'.
    - More details, to see document of factor-analyzer.
- `--FA_method`: To choose method of factor analysis; default is 'minres'.
    - You can choose 'minres', 'ml' or 'principal'.
    - More details, to see document of factor-analyzer.
- `--FA_n_factor`: Number of factors in factor analysis; default is 3.
    - If $kaiser_guttman is activate, this argument is ignored.
- `--kaiser_guttman`: If activate, number of factors is decided based on kaiser-guttman method automatically.
- `--fig_ext`: To choose output figures' extention; default is 'png'.
    - You can choose 'png', 'eps' or 'pdf'.
- `--fig_dpi`: Dpi of output figures; defautl is 300.
- `--heatmap_figsize_x`: X-size of output heatmap that is shows factor loadings; default is 5.
- `--heatmap_figsize_y`: Y-size of output heatmap that is shows factor loadings; default is 6.
- `--heatmap_sort_factor`: To choose sort criteria of heatmap; default is 1 (i.e., Factor 1).
- `--heatmatp_cmap`: To choose heatmaps' color map; default is None (i.e., to be decided automatically).
    - More details, to set dotument of matpltolib.

## Output Files

- eigen_value.png
    - This figure shows eigen values of factor analysis.
    - If activate $kaiser_guttman, this figure is saved.
- factor_loadings.tsv
    - In this file, factor loadings are saved.
- factor_loadings_heatmap.png
    - This figure shows factor loadings as heatmap.
- factor_score.tsv
    - In this file, factor scores are saved.
- rslt.txt
    - Log

Figures' extention is decided by $fig_ext.