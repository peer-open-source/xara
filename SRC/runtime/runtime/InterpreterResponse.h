class Response;

class InterpreterResponse {
public:
  enum class Type {
    Vector,
    Double
  };
  InterpreterResponse(Response* response, Type type, int size, double* data)
    : response(response), type(type), size(size), data(data) {}

  ~InterpreterResponse() {
    delete response;
  }

public:
  int size;
  Type type;
  double* data;
  Response* response;
};