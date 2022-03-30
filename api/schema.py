import federation

from .query import Query, Types
from .mutation import Mutation


SCHEMA = federation.Schema(
    query=Query,
    #mutation=Mutation,
    types=Types
)
