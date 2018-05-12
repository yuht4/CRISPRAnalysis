#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MySQL Database connection, Oracle DB connection

use jgi_connect_db(db_name):
rqc-dev
rqc-prod
sdm
dw
img
seq-location-repos

"""

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## libraries to use


import MySQLdb
import cx_Oracle
import base64 # encrypt/decrypt mysql user's pwd
import json
import os
import time # for timestamp
import getpass # logging users connecting to rqc


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function definitions

'''
create database connection
- uses defaults if nothing passed

To be deprecated: BF 2017-11-27
* use jgi_connect_db instead
'''
def db_connect(db_server = None, db_name = None, db_user = None, db_pwd = None):

    if not db_name:
        db_name = "brycef_rqc"

    if not db_user:
        db_user = "rqc_read"

    if not db_pwd:
        db_pwd = "d29iZWdvbg=="

    if not db_server:
        db_server = "spin"

    if db_server == "dev":
        db_server = "spin"
        #db_server = "draw.jgi-psf.org"
        #db_name = "rqc_dump"

    if db_server == "spin":
        db_server = "db.rqcdev.prod-cattle.stable.spin.nersc.org"
        db_name = "rqc"
        db_pwd = "c25nUjYhdm1CenFjS3k2JQ=="

    if db_server in ("gpdb09.nersc.gov", "production"):
        db_server = "nerscdb04.nersc.gov"
        #db_server = "scidb1.nersc.gov"
        #db_server = "gpdb09.nersc.gov"


    db = None
    db_pwd = base64.decodestring(db_pwd)


    conn_cnt = 0 # how many connection attempty=s
    done = 0
    while done == 0:

        conn_cnt += 1

        try:
            #print "%s: Connecting to %s@%s" % (conn_cnt, db_name, db_server)
            db = MySQLdb.connect(host = db_server, user = db_user, passwd = db_pwd, db = db_name, charset = "utf8", use_unicode = True)
            db.autocommit(True) # so we don't need to commit each time

            done = 1

        except MySQLdb.Error, e:
            print "Error: DB Connection Error %d: %s" % (e.args[0], e.args[1])

        if conn_cnt >= 3:
            done = 1

        # try sleeping 5 seconds and retry to connect
        if done == 0:
            pass



    return db




'''
Open db connection - all read only!
mysql:
rqc dev
rqc prod
sdm (maybe)

oracle:
data warehouse
img

jgi_connect_db("rqc_dev")


# m = mode: q = quiet - don't print anything (default), otherwise print stuff

'''
def jgi_connect_db(db_name, m = "q"):

    db = None
    do_log = True # write to a log file

    if db_name:
        db_name = db_name.lower().strip().replace("_", "-") # rqc_dev -> rqc-dev

        # default engine = mysql

        # file with connection info, not in repository
        j_file = "/global/dna/projectdirs/PI/rqc/prod/versions/jgi_connect_db.json"

        j = {}
        if os.path.isfile(j_file):
            fh = open(j_file, "r")
            j = json.load(fh)
            fh.close()
        else:
            print "Error: cannot find %s" % (os.path.basename(j_file))
            return db

        # aliases
        db_alias = { "prod" : "rqc-prod", "dev" : "rqc-dev", "its" : "dw", "rqc" : "rqc-prod" }
        if db_name in db_alias:
            db_name = db_alias[db_name]




        if db_name in j:

            # log to a file
            if do_log:
                timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
                
                log_file = "/global/projectb/scratch/qc_user/rqc/prod/logs/db/jgi_connect_db.txt" 
                log_buffer = "%s|%s|%s" % (timestamp, getpass.getuser(), db_name)
                
                if not os.path.isfile(log_file):
                    fh = open(log_file, "w")
                    fh.write("#Timestamp|User|Db Connection\n")
                    fh.close()
                    os.chmod(log_file, 0777)
                    
                fh = open(log_file, "a")
                fh.write(log_buffer + "\n")
                fh.close()


            engine = "mysql"
            if 'engine' in j[db_name]:
                engine = str(j[db_name]['engine']).lower()

            if engine == "mysql":
                db_pwd = j[db_name]['pwd']
                try:
                    db = MySQLdb.connect(host = j[db_name]['host'], user = j[db_name]['user'], passwd = db_pwd, db = j[db_name]['db'], charset = "utf8", use_unicode = True)
                    db.autocommit(True) # so we don't need to commit each time
                    if m != "q":
                        print "Connected to MySQL: %s, %s@%s" % (db_name, j[db_name]['user'], j[db_name]['host'])
                except:
                    if m != "q":
                        print "Error: connection to %s failed!" % (db_name)

                    # write to a log somewhere?
                    #sys.exit(4)
            elif engine == "oracle":
                conn_str = "%s/%s@%s" % (j[db_name]['user'], j[db_name]['pwd'], j[db_name]['db'])

                db = cx_Oracle.connect(conn_str)
                if not db:
                    if m != "q":
                        print "Error: connection to %s failed!" % (db_name)
                else:
                    if m != "q":
                        print "Connected to Oracle: %s" % db_name



            else:
                if m != "q":
                    print "Error: unrecognized db engine '%s'" % str(j[db_name]['engine'])

        else:
            if m != "q":
                print "Error: db name '%s' unrecognized" % db_name

    else:
        if m != "q":
            print "Error: no db name"

    return db



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Main Program


if __name__ == "__main__":

    # unit tests
    db = jgi_connect_db("dev", "d")
    if db:
        sth = db.cursor(MySQLdb.cursors.DictCursor)
        sql = "select seq_unit_name, library_name from seq_units where sow_item_id = %s"
        sth.execute(sql, (239289)) # BZZHY, 2 seq units
        row_cnt = int(sth.rowcount)
        for _ in range(row_cnt):
            rs = sth.fetchone()
            print "* %s = %s" % (rs['seq_unit_name'], rs['library_name'])

    else:
        print "No db connection to dev!"

    #db = jgi_connect_db("img", "d")
    db = jgi_connect_db("dw", "d")
    if db:
        sth = db.cursor()
        sql = "select rqc_seq_unit_name, lib_name, sow_status from all_inclusive_report where sow_item_id = :sid"
        sth.execute(sql, sid = 239289)
        for rs in sth:
            print "* %s = %s, %s" % (rs[0], rs[1], rs[2])
    else:
        print "No db connection to dw!"

    exit(0)

