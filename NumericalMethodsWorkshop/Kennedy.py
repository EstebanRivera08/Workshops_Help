

def StringChallenge(st) :
    st = st.replace("%"," ")
    st = st.replace("#"," ")
    st = st.replace("-"," ")
    st = st.replace("+"," ")
    st = st.replace("*"," ")

    lst = st.split()
    st1 = ''
    for i in lst[1:] :
        st1 = st1 + i[0].upper() + i[1:].lower()
    
    result = lst[0].lower() + st1
    print(result)


StringChallenge("BOB loves codifing")
StringChallenge("KENNEDY boBa hIjUePuTa")
StringChallenge("KENNEDY+ho%jij-boBa*hIjUePuTa")

st = "KENNEDY+ho%jij-boBa*hIjUePuTa"
print(st[0])
print(range(st))
