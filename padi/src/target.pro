--  
-- Protocol for Application Database Interface (PADI)
--
-- Copyright 1995, 1996  Bank of Canada.
--
-- The user of this software has the right to use, reproduce and distribute it.
-- Bank of Canada makes no warranties with respect to the software or its 
-- fitness for any particular purpose. The software is distributed by the Bank
-- of Canada solely on an "as is" basis. By using the software, user agrees to 
-- accept the entire risk of using this software.
--

procedure d_target
argument %name
--
-- create test series
--
frequency daily
date 1jan91 to *
case 1 to *
series id(%name + ".c1.daily"):precision by date = 1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1
frequency business
series id(%name + ".c2.business"):precision by date = 10,20,30,40,50,60,70,80,90,80,70,60,50,40,30,20,10
frequency weekly(friday)
series id(%name + ".c3.weekly"):precision by date = 100,200,300,400,500,600,700,800,900,800,700,600,500,400,300,200,100
frequency monthly
series id(%name + ".c3.monthly"):precision by date = 1000,2000,3000,4000,5000,6000,7000,8000,9000,8000,7000,6000,5000,4000,3000,2000,1000
frequency quarterly(december)
series id(%name + ".c4.quarterly"):precision by date = 10000,20000,30000,40000,50000,60000,70000,80000,90000,80000,70000,60000,50000,40000,30000,20000,10000
frequency annual(december)
series id(%name + ".c5.annual"):precision by date = 100000,200000,300000,400000,300000,200000,100000
series id(%name + ".c0.case"):precision by case = 1.1,2.2,3.3,4.4,3.3,2.2,1.1

--
-- create test formula
--
formula id(%name + ".f4.quarterly") = (etsusa.c4.quarterly/10000 + etscdnfor.c4.quarterly)/10000

--
-- create "static var" for test function
--
scalar id(%name + ".f5.count"):numeric = 0

end procedure

