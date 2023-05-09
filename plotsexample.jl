
using Plots
gr()

function main()
    for i in 1:30
        x = randn(10)
        y = randn(10)
        p = scatter(x, y, showaxis=false, legend=false, grid=false, reuse=true, show=true)
        display(p)
        sleep(0.1)
    end
end

main()