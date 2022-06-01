
# client side
import socket
import select
import sys
import datetime

server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
if len(sys.argv) < 3:
    print(
        f"Usage: python3 {__file__} server_IP_address server_port_number, [UserName]")
    exit()

IP_address = str(sys.argv[1])
Port = int(sys.argv[2])
if len(sys.argv) >= 4:
    UserName = "_".join(sys.argv[3:])
else:
    UserName = "Unknown"
server.connect((IP_address, Port))
server.send(UserName.encode())


while True:
    try:
        sockets_list = [sys.stdin, server]
        read_sockets, write_socket, error_socket = select.select(
            sockets_list, [], [])

        for socks in read_sockets:
            if socks == server:
                try:
                    message = socks.recv(2048)
                    print(message.decode())
                except Exception:
                    print("error")
            else:
                message = sys.stdin.readline()
                server.send(message.encode())
    except Exception:
        break

server.close()
