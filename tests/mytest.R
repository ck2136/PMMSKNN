app <- ShinyDriver$new("../")
app$snapshotInit("mytest")

# Input 'table1_rows_current' was set, but doesn't have an input binding.
# Input 'table1_rows_all' was set, but doesn't have an input binding.
app$setInputs(patient_id = 2)
app$setInputs(age_at_s = 74)
app$setInputs(sex1 = "Male")
app$setInputs(submit = "click")
app$setInputs(submit = "click")
# Input 'table1_rows_current' was set, but doesn't have an input binding.
# Input 'table1_rows_current' was set, but doesn't have an input binding.
app$setInputs(delprev = "click")
# Input 'table1_rows_current' was set, but doesn't have an input binding.
# Input 'table1_rows_all' was set, but doesn't have an input binding.
# Input 'table1_rows_current' was set, but doesn't have an input binding.
app$snapshot()
