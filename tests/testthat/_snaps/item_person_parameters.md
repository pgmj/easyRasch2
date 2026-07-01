# output = 'file' errors early when filename is missing

    Code
      RMitemParameters(make_polytomous(), output = "file")
    Condition
      Error:
      ! `filename` must be a single non-empty string when `output = "file"`.

