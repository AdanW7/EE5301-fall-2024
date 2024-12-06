# PA2 -- Adan Wodzinski

	To compile the program simply run 'make' 

	This will compile all of the relevant CPP files
	which can me ran by typing ./placer argument
	where argument is the name of the files we want to examine.
	IE ./placer toy01
	this will parse the neccessary argument files, determine all of the 
	relevant data (Q, X, Y, dx, dy, wire length) and create our CSV's for 
	before and after cell spreading 
	
	This will also create a virtual environment folder for python
	and install all the relevant libraries through pip install.
	to create a images for argument files type into your terminal 
	make run file=argument , where argument is the same as used in placer.
	this will create two matlab plots using the csv data generated from placer.
	
	To simplify the entire process I created a basic script that will allow you to type 
	"./PA2_Tester argument" Ie "./PA2_Tester toy01" and it will run both placer and the python script 
	requires you to have type "make" atleast one time