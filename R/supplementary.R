# Perform do.call on a function `what` but only passing the arguments from
# `args` that match the arguments of `what`. This avoids errors/warnings from
# passing more arguments than necessary
do.call.matched <-
  function(
    what,
    args,
    ...
  ){
    matched_args <-
      intersect(
        x = formalArgs(what),
        y = names(args)
      )

    do.call(what = what, args = args[matched_args], ...)
  }
