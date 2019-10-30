pkgconfig::set_config("drake::strings_in_dots" = "literals")

my_plan <- drake_plan(
  raw_data = GetData(),
  aggregate_data = ResolveNames(FixCommonNames(AggregateData(raw_data))),
  aggregate_csv = write.csv(aggregate_data, file=file_out("data/aggregate_data.csv")),
  all_trees = GetTrees(aggregate_data),
  bullseye_plot = PlotTreeWithTraits(SanitizeTree(all_trees$otol), SanitizeData(aggregate_data)),
  heatmap_plot = PlotPhyloHeatmap(SanitizeTree(all_trees$otol), SanitizeData(aggregate_data)),
  individual_plot = PlotIndividualTraits(SanitizeTree(all_trees$otol), SanitizeData(aggregate_data)),
  corhmm_results = RunCorHMM(SanitizeTree(all_trees$otol), SanitizeData(aggregate_data)),
  save_corhmm = save(corhmm_results, file=file_out("data/corhmm.rda")),
  corhmm_plots = PlotCorHMM(corhmm_results)
)
