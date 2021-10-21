mp.v.modify <-
function (model)
{
    which.term <- names(model$G.param)[grep("X\\.", names(model$G.param))]
    if (length(which.term) > 0) {
        model$G.param <- mpwgaim:::mp.v.init(which.term, model$G.param)
    }
    which.term <- grep("mbf(\"ints\")", names(model$G.param))
    if (length(which.term) > 0) {
        model$G.param <- mpwgaim:::mp.v.init("mbf(\"ints\")", model$G.param)
    }
    model
}
