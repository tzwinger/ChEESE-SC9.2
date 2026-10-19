-- Material parameters for ice

yearinsec=365.25*24*60*60
rhoi=917.0
rhow=1025.0
gravity=9.81
scaledgravity=gravity*(yearinsec^2.0)
scaledrhoi=rhoi/(1.0e06*yearinsec^2.0)
Eyoung=9.0e09
constvisc=1.0e11

-- Elastic/film parameters of the case

E=9.0e9
nu=1.0/2.0
etaw=1.7e-3 -- dynamic viscosity at 0°C
H=200.0 
L=1000.0    -- 1000 or 10000
Ly=500.0    -- or 1000
D=Eyoung*H^3/(12*(1-nu^2))
K=D*24*(16/5)/Ly^4
AC=2*1.7/K
bslp=0.0
Pext=2.0e5
maxy=8.625E06

-- ! * 1.0e-06 * yearinsec

-- Glen's exponent and inverse
n = 3.0
m = 1.0/n

-- rate factors and molar activiation energies
-- after Paterson 2010
A1 = 2.89165e-13
A2 = 2.42736e-02
Q1 = 60.0e3
Q2 = 115.0e3

Tc = -0.1
Arrh = 3.0E-24

offsetpressure=1.0E04
h0=0.001

function tides(t)
   sl = 3.0 + math.sin(4.0*math.pi*t/(24*3600))*(0.5 + 0.5*abs(math.sin(pi*t/(24.0*3600.0))))
   -- print(">>>> Sea level:",sl)
   return sl
end

function seapressure(t,z)
  sl = tides(t)
  if (z > sl) then
       pw = 0.0
  else   
       pw = rhow * gravity * (sl - z)
       -- print(">>>>> pw:",pw)
  end
  return pw
end  

function seapressurestatic(z)
  sl = 5
 if (z > sl) then
       pw = 0.0
  else
       pw = rhow * gravity * (sl - z)
end
return pw
end

function IfThenElse(cond,pos,neg)
   if (cond) then
      return pos
   else
      return neg
   end
end
