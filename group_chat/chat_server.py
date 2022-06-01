
# server side
import socket
import sys
import _thread
import datetime

server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

if len(sys.argv) != 3:
    print(f"Usage: python3 {__file__} [IP address] [port number]")
    exit()

IP_address = str(sys.argv[1])
Port = int(sys.argv[2])
server.bind((IP_address, Port))
server.listen(100)
list_of_clients = []
list_of_userNames = []


def clientthread(conn, addr, UserName):
    with open("chat_history.txt", 'r') as f:
        chat_history = "\n".join(f.readlines())

    # conn.send(chat_history.replace(f"{UserName}<{addr[0]}>:", "You:").encode())
    try:
        print(chat_history.encode())
        conn.send(chat_history.encode())
        conn.send(f"\nWelcome {UserName}!".encode())
    except Exception:
        print('err')

    while True:
        try:
            message = conn.recv(2048)
            if message:
                if message.decode()[0] == "/":
                    msg_list = message.decode().split(" ")
                    if len(msg_list) < 3:
                        conn.send(help().encode())
                    if msg_list[0] in ["/msg", "/pv"]:
                        ign = msg_list[1]
                        message_to_send = ' '.join(msg_list[2:])
                        private_msg(UserName, ign, message_to_send)
                else:
                    # message_to_send = f"{UserName}<{addr[0]}>: {message}"
                    time = datetime.datetime.now().time()
                    message_to_send = f"\n{UserName}:\n\t{message.decode()}"
                    message_to_send += f"\t\t\t{time.hour :02d}:{time.minute:02d}\n"
                    broadcast(message_to_send, conn)
            # else:
            #     remove(conn)
            #     conn.close()
            #     break
        except Exception:
            remove(conn)
            conn.close()
            break


def broadcast(message, connection):
    for client in list_of_clients:
        # if client != connection:
        try:
            client.send(message.encode())
        except Exception:
            client.close()
            remove(client)
    try:
        with open("chat_history.txt", "a") as f:
            print(message, file=f)
    except Exception:
        connection.close()
        remove(connection)
        return


def help():
    return "Usage /msg UserName your massage"


def remove(connection):
    if connection in list_of_clients:
        i = list_of_clients.index(connection)
        UserName = list_of_userNames[i]
        broadcast(f"\t--- {UserName} left the room ---", connection)
        list_of_clients.remove(connection)
        list_of_userNames.remove(UserName)


def private_msg(UserName, ign, message):
    try:
        if ign in list_of_userNames:
            j = list_of_userNames.index(ign)
            conn_to_send = list_of_clients[j]
            message_to_send = f"*****\nPrivate Message from {UserName}:\n{message}\n*******"
            conn_to_send.send(message_to_send.encode())
            msg = f"******\nPrivate Message to {ign}\n{message}\n********"
            i = list_of_userNames.index(UserName)
            conn = list_of_clients[i]
            conn.send(msg.encode())
        else:
            i = list_of_userNames.index(UserName)
            conn = list_of_clients[i]
            conn.send(f"{ign} is not online!".encode())
    except Exception:
        print('error sending private message')


while True:
    try:
        conn, addr = server.accept()

        UserName = conn.recv(1024).decode()
        if UserName in list_of_userNames:
            UserName = f"{UserName}<{addr[0]}>"

        list_of_clients.append(conn)
        list_of_userNames.append(UserName)
        message = f"\t--- {UserName} connected to the room ---"
        broadcast(message, conn)
        _thread.start_new_thread(
            clientthread, (conn, addr, UserName))
    except Exception:
        break

conn.close()
server.close()
