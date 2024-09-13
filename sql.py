import mysql.connector as con

connection = con.connect(
      host='localhost',
      user='root',
      password='root',
      database='checked_functions'
)

mycursor = connection.cursor()


sql2 = "INSERT INTO functions (Phi_x, Phi_y, ExtensionDegree, CheckResult, CheckAmount) VALUES ('x*y', 'y', -1, 'CONTRACTS', 0);"
sql = "SELECT * FROM checked_functions.functions;"
mycursor.execute(sql2)
connection.commit()
result = mycursor.fetchall()
print(result)