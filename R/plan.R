pkgconfig::set_config("drake::strings_in_dots" = "literals")

my_plan <- drake_plan(
  raw_data = GetData(),
  aggregate_data = AggregateData(raw_data)
)
