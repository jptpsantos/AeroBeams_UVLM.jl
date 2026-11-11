using Printf
mkpath(joinpath(pwd(),"tmp","hub_pitch_test"))
path=joinpath(pwd(),"tmp","hub_pitch_test","synthetic_history.txt")
open(path,"w") do io
    println(io,"Synthetic test history")
    println(io,"time_s\tpropeller_1_pitch_deg\tpropeller_1_yaw_deg\tpropeller_1_attachment_wing_twist_deg\ttip_displacement_m\ttip_twist_deg")
    for t in 0.0:0.002:5.0
        a=2π*4.0*(t-3.95)
        envelope=exp(0.05*(t-3.95))
        wing=0.12*envelope*cos(a)
        attach=0.1*envelope*cos(a+deg2rad(10.0))
        relative=2.0*envelope*cos(a+deg2rad(30.0))
        yaw=2.1*envelope*cos(a+deg2rad(120.0))
        @printf(io,"%.6f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n",t,relative,yaw,attach,0.0,wing)
    end
end
println(path)
