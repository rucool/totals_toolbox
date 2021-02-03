function dbConn = mongo_database
    %Connect to MongoDB database
    port = 27017;
    dbname = 'codar';
  
    host = '127.0.0.1';
    username = 'admin';
    password = 'password';
    dbConn = mongo(host, port, dbname, "UserName", username, "Password", password);
end