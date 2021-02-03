function dbConn = mysql_database
    %Connect to MySQL cool_ais database
    host = 'localhost';
    user = 'admin';
    password = 'root';
    dbName = 'coolops';
    jdbcString = sprintf('jdbc:mysql://%s/%s', host, dbName);
    jdbcDriver = 'com.mysql.jdbc.Driver';
    % Make the connection
    dbConn = database(dbName, user , password, jdbcDriver, jdbcString);
end