exp=math.exp
abs=math.abs
sin=math.sin
pi=math.pi
sqrt=math.sqrt
T0=0.0
-- remember LUA is case sensitive! MAKE SURE TEMPERATURE BOTH IN CELSIUS HERE AND SIF?? CHECK??
-- NOTE TESTED OTHER TEMPERATURE PROFILE HERE!!
function temperature(d)
 -- print("temperature: d = ", d)
 z = -d   -- some comment; z is minus depth
--if (1<0) then 
 if (z<-260) then
  return T0
 elseif (z<-210) then
  return linint(-260.0,T0,-210.0,-3.5,z)
 elseif (z<-150) then
  return linint(-210.0,-3.5,-150.0,-6.0,z)
 elseif (z<-90) then
  return linint(-150.0,-6.0,-90.0,-6.0,z)
 elseif (z<-25) then
  return linint(-90.0,-6.0,-25.0,-2.5,z)
 elseif (z<0) then
  return linint(-25.0,-2.5,0.0,-10.0,z)
 else
  return -10.0
 end
--else 
-- return -25.0 + 20.0*d / 1500.0
--end 
end

function linint(z1,t1,z2,t2,z)
  return  t2*(z-z1)/(z2-z1) + t1*(z2-z)/(z2-z1) 
end


