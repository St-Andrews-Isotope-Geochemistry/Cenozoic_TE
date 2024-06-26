from pathlib import Path
import os
import pandas as pd
import numpy as np
import sqlalchemy as sqla
from sqlalchemy.exc import SQLAlchemyError
#from sqlalchemy.orm import sessionmaker


## functions

def query_to_dataframe(sql_query):
    """
    Executes an SQL query and returns the result as a pandas DataFrame.
    
    Parameters:
    - database_url (str): The database connection URL.
    - sql_query (str): The SQL query to execute.
    
    Returns:
    - pd.DataFrame: The query result as a DataFrame, or None if an error occurs.
    """
    # Create the SQLAlchemy engine
    #engine = create_engine(database_url)

    try:
        # Execute the query and fetch the results
        with engine.connect() as conn:
            result = conn.execute(sqla.text(sql_query))

            # Convert the query result to a pandas DataFrame
            df = pd.DataFrame(result.fetchall(), columns=result.keys())

        return df
    except SQLAlchemyError as e:
        # Print the error message
        print(f"An error occurred: {e}")
        return None


def execute_query(sql_query):
    """
    Executes an SQL query and returns the result as a pandas DataFrame.
    
    Parameters:
    - database_url (str): The database connection URL.
    - sql_query (str): The SQL query to execute.
    
    Returns:
    - pd.DataFrame: The query result as a DataFrame, or None if an error occurs.
    """
    # Create the SQLAlchemy engine
    #engine = create_engine(database_url)

    try:
        # Execute the query and fetch the results
        with engine.connect() as conn:
            result = conn.execute(sqla.text(sql_query))

        return result
    except SQLAlchemyError as e:
        # Print the error message
        print(f"An error occurred: {e}")
        return None



#define paths
data_path=Path(os.getcwd())/"data"


foram_df=pd.read_csv(data_path/"foram_dataframe.csv", index_col=0)

DATABASE_URL="sqlite:///data/foram_database.db"
engine = sqla.create_engine(DATABASE_URL)
#Session = sessionmaker(bind=engine)
#session = Session()
conn = engine.connect()

foram_df.to_sql('master', con=engine, if_exists='replace', index=False)



#Create core table from master
execute_query("""CREATE TABLE cores 
                AS SELECT DISTINCT core AS core_id, lat AS latitude, lon AS longitude, ocean_basin
                FROM master""")

TE_list=['Li7', 'B11', 'Na23', 'Mg24', 'Mg25', 'Mg26', 'Al27', 
         'P31', 'S32', 'K39', 'Ni58', 'Co59', 'Cu63', 'Zn64', 'Rb85',
         'Sr88', 'Mo92', 'Cd111', 'Ba138', 'Nd146', 'U238']

TE_list_joined=', '.join(TE_list)
condition = " IS NOT NULL OR ".join(TE_list)+" IS NOT NULL"


#setup trace_elements table with primary key
execute_query(f"""CREATE TABLE trace_elements (
                te_id INTEGER PRIMARY KEY, 
                {TE_list_joined})""")

#copy over the data from master
execute_query(f"""INSERT INTO trace_elements ({TE_list_joined})
                SELECT {TE_list_joined}
                FROM master
                WHERE {condition}""")

#retrieve trace_elements table
df = query_to_dataframe("SELECT * FROM trace_elements")
