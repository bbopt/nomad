function stop = cbFun(count,f,x)

stop = false;
if (count > 10 && f < 1.0 )
    stop = true;
end
return
