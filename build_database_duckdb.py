from pathlib import Path
import os
import pandas as pd
import numpy as np
import duckdb 
#import sqlalchemy as sqla
#from sqlalchemy.exc import SQLAlchemyError
#from sqlalchemy.orm import sessionmaker



#Add column with FK to table UNFINISHED
def add_fk(table, column, fk_table, fk_column):
    
    #find the constraints of the current table
    constraints_table=duckdb.sql(f"""PRAGMA table_info({table})""").df()
    constraints_table['fk']=False
    constraints=duckdb.sql(f"""SELECT kcu1.table_name AS t1, kcu1.column_name AS col1, kcu2.table_name AS t2, kcu2.column_name AS col2
    FROM information_schema.key_column_usage AS kcu1
    INNER JOIN information_schema.referential_constraints AS rc
    ON kcu1.constraint_name = rc.constraint_name
    INNER JOIN information_schema.key_column_usage AS kcu2
    ON rc.unique_constraint_name = kcu2.constraint_name
    WHERE kcu1.table_name = {table}""").df()
    for i, row in constraints.iterrows():
        idx=constraints_table['name']==row['col1']
        constraints_table.loc[idx, 'fk']=True
        constraints_table.loc[idx, 'fk_tab']=row['t2']
        constraints_table.loc[idx, 'fk_col']=row['c2']
    
    #create the new table with the new column
    duckdb.sql(f"""CREATE TABLE dummy_tab (
    core_id VARCHAR(50) PRIMARY KEY, 
    latitude REAL, 
    longitude REAL, 
    ocean_basin VARCHAR(50))""")
    
    
    
    
    query=f"""ALTER TABLE {table}
    ADD FOREIGN KEY ({column})
    REFERENCES {fk_table} ({fk_column})"""
    duckdb.sql(query)



## paths
#define paths
data_path=Path(os.getcwd())/"data"
database_path=Path(os.getcwd())/"foram_database"

master_file=data_path/"foram_dataframe_240625.csv"

foram_df=pd.read_csv(master_file, index_col=0)
database_path

# Connect to an in-memory DuckDB database
conn = duckdb.connect(database=str(database_path/"foram_database.db"))

# Load master data into a table
query = f""" CREATE TABLE master AS
SELECT * 
FROM read_csv('{master_file}')
"""
duckdb.sql(query)

## Core table

#Setup core table from master
duckdb.sql("""CREATE TABLE cores (
    core_id VARCHAR(50) PRIMARY KEY, 
    latitude REAL, 
    longitude REAL, 
    ocean_basin VARCHAR(50))""")

#insert data into core table
duckdb.sql("""INSERT INTO cores (core_id, latitude, longitude, ocean_basin)
                SELECT DISTINCT core AS core_id, lat AS latitude, lon AS longitude, ocean_basin
                FROM master""")


#Retrieve core table
core_table=duckdb.sql("SELECT * FROM cores").df()

## Trace elements table

TE_list=['Li7', 'B11', 'Na23', 'Mg24', 'Mg25', 'Mg26', 'Al27', 
         'P31', 'S32', 'K39', 'Ni58', 'Co59', 'Cu63', 'Zn64', 'Rb85',
         'Sr88', 'Mo92', 'Cd111', 'Ba138', 'Nd146', 'U238']

TE_list_joined=', '.join(TE_list)
condition = " IS NOT NULL OR ".join(TE_list)+" IS NOT NULL"


#setup trace_elements table with PRIMARY KEY
duckdb.sql(f"""CREATE TABLE trace_elements (
    te_id TEXT PRIMARY KEY,
    core_id VARCHAR(100) REFERENCES cores(core_id),
    species VARCHAR(50),
    {' REAL, '.join(TE_list)} REAL)""")


# insert into trace_elements table
duckdb.sql(f"""
INSERT INTO trace_elements (te_id, core_id, species, {TE_list_joined})
SELECT 
    sample_mix_id || '_' || species AS te_id, 
    core AS core_id, 
    species,
    {', '.join([f'AVG({col}) AS {col}' for col in TE_list])}
FROM master
WHERE {condition}
GROUP BY 
    sample_mix_id,
    species, 
    core
""")


#retrieve trace_elements table
te_df=duckdb.sql("SELECT * FROM trace_elements").df()


## d11B table

#setup d11B table
duckdb.sql("""CREATE SEQUENCE d11B_id_sequence START 1;
           CREATE TABLE boron_isotopes (
               d11B_id INTEGER PRIMARY KEY DEFAULT nextval('d11B_id_sequence'),
               core_id VARCHAR(50) REFERENCES cores(core_id),
               species VARCHAR(50),
               d11B REAL, 
               d11B_2se REAL)""")

#insert data into d11B table
duckdb.sql("""INSERT INTO boron_isotopes (core_id, species, d11B, d11B_2se)
                SELECT core AS core_id, species, d11B, d11B_2SE AS d11B_2se
                FROM master
                WHERE d11B IS NOT NULL OR d11B_2se IS NOT NULL""")


## sample_info_table (age models)

#setup age_models table
duckdb.sql("""CREATE SEQUENCE age_model_id_sequence START 1;
           CREATE TABLE age_models (
               age_model_id INTEGER PRIMARY KEY DEFAULT nextval('age_model_id_sequence'),
               core_id VARCHAR(50) REFERENCES cores(core_id),
               age_Ma REAL,
               """)


