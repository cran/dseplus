--  
-- Protocol for Application -- Database Interface (PADI)
--
-- Copyright 1995, 1996  Bank of Canada.
--
-- The user of this software has the right to use, reproduce and distribute it.
-- Bank of Canada makes no warranties with respect to the software or its 
-- fitness for any particular purpose. The software is distributed by the Bank
-- of Canada solely on an "as is" basis. By using the software, user agrees to 
-- accept the entire risk of using this software.
--

PROCEDURE D_OBJECT
---------------------------------------------
-- create object "tables" in current database
---------------------------------------------

scalar obj_dbdir:string 
scalar obj_sysdbpath:string
description(obj_dbdir) = "Database Directory"
description(obj_sysdbpath) = "System Database Path"

series obj.name:string by case
series obj.desc:string by case
series obj.dbname:string by case
series obj.access:string by case
series obj.class:string by case
series obj.freq:string by case
series obj.type:string by case
series obj.fame_name:string by case
series obj.fame_class:string by case
series obj.fame_ldpath:string by case
description(obj.name) = "Object Name"
description(obj.desc) = "Object Description"
description(obj.dbname) = "Object Database Name"
description(obj.access) = "Object Access Mode (READ, UPDATE, SHARED)"
description(obj.class) = "Object Class (SERIES or SCALAR)"
description(obj.freq) = "Object Frequency"
description(obj.type) = "Object Data Type"
description(obj.fame_name) = "Fame Object Name"
description(obj.fame_class) = "Fame Object Class (SERIES, SCALAR, FORMULA, FUNCTION)"
series obj_link.dbname:string by case
series obj_link.access:string by case
description(obj_link.dbname) = "Linked Database Name"
description(obj_link.access) = "Linked Database Access Mode (READ, UPDATE, SHARED)"

END PROCEDURE

PROCEDURE S_DBDIR
ARGUMENT %dbdir
----------------------------------------
-- set the database directory
----------------------------------------
if not exists(obj.name)
  d_object
end if

set obj_dbdir = %dbdir

END PROCEDURE

PROCEDURE S_SYSDBPATH
ARGUMENT %sysdbpath
----------------------------------------
-- set the database directory
----------------------------------------
if not exists(obj.name)
  d_object
end if

set obj_sysdbpath = %sysdbpath

END PROCEDURE

PROCEDURE ADD_OBJ
ARGUMENT %name, %desc, %dbname, %access, %class, %freq, %type, %fame_name, %fame_class, %fame_ldpath
----------------------------------------
-- add an object to the current database
----------------------------------------
if not exists(obj.name)
  d_object
end if

local'%i = LASTVALUE(obj.name) + 1
if missing(%i)
  set %i = 1
end if

set obj.name[%i] = %name
set obj.desc[%i] = %desc
set obj.dbname[%i] = %dbname
set obj.access[%i] = %access
set obj.class[%i] = %class
set obj.freq[%i] = %freq
set obj.type[%i] = %type
set obj.fame_name[%i] = %fame_name
set obj.fame_class[%i] = %fame_class
set obj.fame_ldpath[%i] = %fame_ldpath

scalar id(%name + ".obj"):numeric = %i
description(id(%name + ".obj")) = "Index into Object Table"

END PROCEDURE

PROCEDURE ADD_DB
ARGUMENT %dbname, %access
--------------------------------------------
-- link a database into the object database
--------------------------------------------
if not exists(obj.name)
  d_object
end if

local'%found = 0
local'%i = LASTVALUE(obj_link.dbname) + 1
if missing(%i)
  set %i = 1
else
-- check to see if database has already been added
  loop for j = 1 to (%i - 1)
     if obj_link.dbname[j] eq %dbname
        set %found = 1
	escape
     endif
  endloop
end if

if %found eq 0
  set obj_link.dbname[%i] = %dbname
  set obj_link.access[%i] = %access
endif

END PROCEDURE

PROCEDURE EX_OBJ
ARGUMENT %object_name
--------------------------------------------
-- exclude an object from a database link
--------------------------------------------
if not exists(obj.name)
  d_object
end if

scalar id(%object_name + ".obj"):numeric = 0
description(id(%object_name + ".obj")) = "Null index into object table for excluded object name"

END PROCEDURE
