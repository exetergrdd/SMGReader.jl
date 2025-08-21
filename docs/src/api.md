# API

- [External API](#External-API)
- [Internal API](#Internal-API)
- [Base extensions](#Base-Extensions)

## External API
```@autodocs
Modules = [SMGReader]
Private = false
```

## Internal API
```@autodocs
Modules = [SMGReader]
Private = true
Public = false
Filter = obj -> parentmodule(obj) != Base
```


## Base Extensions
```@autodocs
Modules = [SMGReader]
Private = true
Public = false
Filter = obj -> parentmodule(obj) == Base
```