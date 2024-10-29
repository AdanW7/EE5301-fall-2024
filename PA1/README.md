# EE5301-ProjA1


## `Makefile`
- The `Makefile` compiles the project into an executable:
  - Typing `make` will produce an executable called `sta`.
  - It includes a `test` target, allowing us to run tests by typing `make test`, which will test any circuit and provide relevant gate information as arguments.


## `updated_tester.sh` 
-   A bash script based on previous project solution
    - designed to test the project and 
    - after we run a test copy the data in "[(ckt_traversal.txt)]" 
    - to the corresponding file in the <em> <ins> results </ins></em> directory
    
### `Main.cpp`
- The main entry point that coordinates the execution of helper classes and files.

###  `Circuit.cpp`
- Manages the circuit and stores pointers to `CircuitNode` objects, enabling us to print detailed circuit node information.


### `CircuitNode.cpp`
- Defines setter and getter functions for manipulating and accessing the properties of a circuit node, 
- including its ID, gate type, input/output pads, fan-in/out list, gate information, etc...


### `GateDatabase.cpp`
- Parses the gate library (database) file to retrieve gate-specific information for use in circuit analysis.


## <em> <ins> test </ins></em> directory 

- Contains provided circuit files and the gate library (database) file


## <em> <ins> results</ins></em> directory 
- This directory is created by the `updated_tester.sh` script and stores the results of the tests 