test_that("events_to_net", {
    events_list <- list(
        i = c(1, 2, 3, 1),
        j = c(2, 3, 1, 4),
        t = c(5, 10, 15, 20)
    )
    net <- events_to_net(events_list)
    expect_equal(network::network.size(net),
                 4,
                 tolerance = 0)

})
test_that("events_to_bipartite_net", {
    events_list <- list(
        i = c("A", "B", "C", "D","D"),
        j = c(2, 3, 1, 2, 1),
        t = c(5, 10, 15, 20, 21)
    )
    net <- events_to_bipartite_net(events_list)
    expect_equal(network::network.size(net),
                 7,
                 tolerance = 0)

})
