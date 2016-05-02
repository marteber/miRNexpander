# INSTALLATION
1. Requirements  
  * Python 2.7+ with modules:
     - biopython
     - matplotlib
     - mysqldb
     - networkx 
     - openpyxl 1.9+
     - webcolors
  * MySQL 5.5+
  * *optional:* libSBML 5.10+ built with python API (requirements: swig, libxml2)
   
    Suggested command(s) for installation on Ubuntu 14.04
  ```shell
  $ sudo apt-get install build-essential mysql-server python-biopython \
  python-webcolors python-networkx python-mysqldb python-matplotlib \
  libxml2 libxml2-dev python-dev python-pip zlib1g zlib1g-dev bzip2 libbz2-dev
  $ sudo pip install python-libsbml
  $ sudo pip install openpyxl
  ```

1. Suggestions
  * phpmyadmin (to more easily inspect database)
  * *optional:* Grant FILE permissions to database user  
    **Note: Mysql clients connecting to the database server need to enable infile_local to make use of the above!**

    Suggested command(s) for installation on Ubuntu 14.04
  ```shell
  $ sudo apt-get install phpmyadmin
  $ mysql -u root -p
  > GRANT FILE ON *.* TO <user>@'%';
  > quit
  ```
